extern crate time;
use clap::Parser;
use dbghap::file_reader;
use dbghap::dbg;
use dbghap::parse_cmd_line;
use dbghap::utils_frags;
use std::fs;
use std::path::Path;
use std::time::Instant;
use dbghap::types_structs::*;

//This makes statically compiled musl library
//much much faster. Set to default for x86 systems...
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[allow(deprecated)]
fn main() {
    #![allow(warnings)]

    let options = parse_cmd_line::Options::parse();
    if options.trace {
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Trace)
            .init()
            .unwrap();
    } else if options.debug{
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Debug)
            .init()
            .unwrap();
    }
    else{
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Info)
            .init()
            .unwrap();
    }

    let contig_out_dir = format!("{}", options.output_dir);

    if Path::new(&contig_out_dir).exists() && options.overwrite {
        let res = fs::remove_dir_all(&contig_out_dir);
        if res.is_err() {
            log::error!(
                "Could not remove {} successfully. Exiting.",
                &contig_out_dir
            );
            std::process::exit(1);
        }
    }
    else if Path::new(&contig_out_dir).exists() && !options.overwrite {
        log::error!(
            "Output directory {} already exists. Use --overwrite to overwrite.",
            &contig_out_dir
        );
        std::process::exit(1);
    }

    fs::create_dir_all(&contig_out_dir).unwrap();


    let start_t_initial = Instant::now();
    log::info!("Preprocessing VCF/Reference");
    let start_t = Instant::now();
    let contigs_to_phase;
    contigs_to_phase = file_reader::get_contigs_to_phase(&options.bam_file);
    let mut main_bam = file_reader::get_bam_readers(&options);
    log::debug!("Read BAM file successfully.");

    let mut chrom_seqs = None;

    let vcf_profile = file_reader::get_vcf_profile(&options.vcf_file, &contigs_to_phase);
    let snp_to_genome_pos_map = file_reader::get_genotypes_from_vcf_hts(options.vcf_file.clone());
    log::debug!("Read VCF successfully.");
    if options.reference_fasta != "" {
        chrom_seqs = Some(file_reader::get_fasta_seqs(&options.reference_fasta));
        log::debug!("Read reference fasta successfully.");
    }
    log::debug!("Finished preprocessing in {:?}", Instant::now() - start_t);

    let mut warn_first_length = true;
    for contig in contigs_to_phase.iter() {

        if options.sequences_to_phase.is_some() && !options.sequences_to_phase.as_ref().unwrap().contains(&contig)
        {
            continue;
        } else if !vcf_profile.vcf_pos_allele_map.contains_key(contig.as_str())
            || vcf_profile.vcf_pos_allele_map[contig.as_str()].len() < options.snp_count_filter
        {
            if warn_first_length {
                log::warn!(
                    "A contig ({}) is not present or has < {} variants. This warning will not be shown from now on. Make sure to change --snp-count-filter if you want to phase small contigs.",
                    contig,
                    options.snp_count_filter,
                );
            }
            warn_first_length = false;
            continue;
        }

        let start_t = Instant::now();
        //log::info!("-----{}-----", contig);
        log::info!(
            "Phasing contig {}",
            contig
        );
        let mut all_frags;
        let frags_without_snps;
        (all_frags, frags_without_snps) = file_reader::get_frags_from_bamvcf_rewrite(
            &mut main_bam,
            &vcf_profile,
            &options,
            &mut chrom_seqs,
            &contig,
        );
        log::debug!("Number of reads passing filtering: {}", all_frags.len());
        if all_frags.len() == 0 {
            log::debug!("Contig {} has no fragments", contig);
            continue;
        }

        if snp_to_genome_pos_map.contains_key(contig) {
            
            let snp_to_genome_pos: &Vec<usize>;
            snp_to_genome_pos = snp_to_genome_pos_map.get(contig).unwrap();

            all_frags.sort();
            for (i, frag) in all_frags.iter_mut().enumerate() {
                frag.counter_id = i;
            }

            //Get last SNP on the genome covered over all fragments.
            let length_gn = utils_frags::get_length_gn(&all_frags);
            log::debug!("Contig {} has {} SNPs", contig, length_gn);

            let mut final_frags;
            final_frags = all_frags;
            let dbg_frags : Vec<FragDBG> = final_frags.iter().map(|x| dbg::frag_to_dbgfrag(x, &options)).collect();

            let dbg_frags = dbg_frags[0..options.max_frags.min(dbg_frags.len())].to_vec();
            log::debug!(
                "Reading inputs, realigning time taken {:?}",
                Instant::now() - start_t
            );

            dbg::construct_dbg(&dbg_frags, &options, &snp_to_genome_pos, &contig);
        }
    }
    log::info!("Total time taken is {:?}", Instant::now() - start_t_initial);
}
