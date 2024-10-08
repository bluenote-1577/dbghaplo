use crate::constants;
use crate::types_structs::FragDBG;
use crate::types_structs::Frag;
use crate::types_structs::{Genotype, GenotypeCount, Haplotype, SnpPosition};
use crate::types_structs::{HapBlock, GAP_CHAR};
use fxhash::{FxHashMap, FxHashSet};
use itertools::Itertools; // 0.8.2
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use statrs::distribution::ChiSquared;
use statrs::distribution::ContinuousCDF;
use std::sync::Mutex;
use std::time::Instant;
use std::mem;

// Get the number # of different bases between the
// two fragments
pub fn distance(r1: &Frag, r2: &Frag) -> (i32, i32) {
    let mut diff = 0;
    let mut same = 0;

    for pos in r1.positions.intersection(&r2.positions) {
        if r1.seq_dict.get(pos) == r2.seq_dict.get(pos) {
            same += (phred_scale(r1, pos) * phred_scale(r2, pos)).round() as i32;
        } else {
            diff += (phred_scale(r1, pos) * phred_scale(r2, pos)).round() as i32;
        }
    }

    (same, diff)
}

pub fn distance_read_haplo_epsilon_empty(r: &Frag, hap: &Haplotype, epsilon: f64) -> (f64, f64) {
    let mut diff = 0.0;
    let mut same = 0.0;
    for pos in r.positions.iter() {
        let mut empty_pos = true;
        if hap.contains_key(pos) {
            for (_key, val) in hap[&pos].iter() {
                if *val != 0. {
                    empty_pos = false;
                    break;
                }
            }
        }
        if empty_pos {
            diff += epsilon;
            continue;
        }

        //Can speed this up by storing the consensus var in Haplotype
        //without need to recompute everytime. 
        let frag_var = r.seq_dict.get(pos).unwrap();
        let consensus_var = hap
            .get(pos)
            .unwrap()
            .iter()
            .max_by_key(|entry| entry.1)
            .unwrap()
            .0;
        if *frag_var == *consensus_var {
            same += phred_scale(r, pos).into_inner();
        } else {
            let frag_var_count = hap.get(pos).unwrap().get(frag_var);
            if let Some(count) = frag_var_count {
                if count == hap.get(pos).unwrap().get(consensus_var).unwrap() {
                    same += phred_scale(r, pos).into_inner();
                    continue;
                }
            }
            diff += phred_scale(r, pos).into_inner();
        }
    }

    (same, diff)
}

pub fn distance_read_haplo(r1: &Frag, hap: &Haplotype) -> (usize, usize) {
    let mut diff = 0.;
    let mut same = 0.;
    for pos in r1.positions.iter() {
        if !hap.contains_key(pos) {
            continue;
        }

        let frag_var = r1.seq_dict.get(pos).unwrap();
        let consensus_var = hap
            .get(pos)
            .unwrap()
            .iter()
            .max_by_key(|entry| entry.1)
            .unwrap()
            .0;
        if *frag_var == *consensus_var {
            same += phred_scale(r1, pos).into_inner();
        } else {
            let frag_var_count = hap.get(pos).unwrap().get(frag_var);
            if let Some(count) = frag_var_count {
                if count == hap.get(pos).unwrap().get(consensus_var).unwrap() {
                    phred_scale(r1, pos).into_inner();
                    continue;
                }
            }
            diff += phred_scale(r1, pos).into_inner();
        }
    }

    (same.round() as usize, diff.round() as usize)
}

//Index each position by the set of fragments which overlap that position
pub fn get_all_overlaps(frags: &Vec<Frag>) -> FxHashMap<SnpPosition, FxHashSet<&Frag>> {
    let mut overlaps = FxHashMap::default();
    for frag in frags.iter() {
        for pos in frag.positions.iter() {
            let pos_set = overlaps.entry(*pos).or_insert(FxHashSet::default());
            pos_set.insert(frag);
        }
    }

    overlaps
}

//Find all distances between fragments.
//Assumes sorted fragments by first position.
pub fn get_all_distances(frags: &Vec<Frag>) -> FxHashMap<&Frag, FxHashMap<&Frag, i32>> {
    let mut pairwise_distances = FxHashMap::default();

    for (i, frag1) in frags.iter().enumerate() {
        let mut from_frag_distance = FxHashMap::default();
        for j in i..frags.len() {
            let frag2 = &frags[j];

            if frag1.last_position < frag2.first_position {
                break;
            }

            from_frag_distance.insert(frag2, distance(frag1, frag2).1);
        }
        pairwise_distances.insert(frag1, from_frag_distance);
    }

    pairwise_distances
}

pub fn check_overlap(r1: &Frag, r2: &Frag) -> bool {
    if r1.last_position < r2.first_position {
        return false;
    }
    if r2.last_position < r1.first_position {
        return false;
    }
    let t: Vec<_> = r1.positions.intersection(&r2.positions).collect();
    if t.len() == 0 {
        return false;
    } else {
        return true;
    }
}

pub fn get_consensus_seq_dict(haplo: &Haplotype) -> FxHashMap<SnpPosition, Genotype> {
    let mut consensus_seq_dict = FxHashMap::default();
    for (pos, geno_dict) in haplo.iter() {
        let consensus_var = geno_dict.iter().max_by_key(|entry| entry.1).unwrap().0;
        consensus_seq_dict.insert(*pos, *consensus_var);
    }
    return consensus_seq_dict;
}

pub fn fragdbg_to_seq_dict(frag_set: &[&FragDBG], use_phred: bool) -> Haplotype {
    let mut hap_map = FxHashMap::default();
    for frag in frag_set.iter() {
        for pos in frag.seq_dict.keys(){
            let var_at_pos = frag.seq_dict.get(pos).unwrap();
            let sites = hap_map.entry(*pos).or_insert(FxHashMap::default());
            let site_counter = sites.entry(*var_at_pos).or_insert(OrderedFloat(0.));
            if use_phred {
                *site_counter += phred_scale_dbg(frag, pos);
            } else {
                *site_counter += 1.;
            }
        }
    }
    return hap_map;
}

pub fn set_to_seq_dict(frag_set: &FxHashSet<&Frag>, use_phred: bool) -> Haplotype {
    let mut hap_map = FxHashMap::default();
    for frag in frag_set.iter() {
        for pos in frag.positions.iter() {
            let var_at_pos = frag.seq_dict.get(pos).unwrap();
            let sites = hap_map.entry(*pos).or_insert(FxHashMap::default());
            let site_counter = sites.entry(*var_at_pos).or_insert(OrderedFloat(0.));
            if use_phred {
                *site_counter += phred_scale(frag, pos);
            } else {
                *site_counter += 1.;
            }
        }
    }
    return hap_map;
}

pub fn hap_block_from_partition(part: &Vec<FxHashSet<&Frag>>, use_qual: bool) -> HapBlock {
    let mut block_vec = Vec::new();
    for set in part.iter() {
        let hap_map = set_to_seq_dict(set, use_qual);
        block_vec.push(hap_map);
    }
    HapBlock { blocks: block_vec }
}

pub fn get_avg_length_dbgf(all_frags: &Vec<FragDBG>, quantile: f64) -> usize{
    let mut length_vec = Vec::new();
    for frag in all_frags.iter() {
        length_vec.push(frag.seq_dict.len() as u32);
    }
    length_vec.sort();
    return length_vec[(length_vec.len() as f64 * quantile) as usize] as usize;
}


pub fn get_avg_length(all_frags: &Vec<Frag>, quantile: f64) -> SnpPosition {
    let mut length_vec = Vec::new();
    for frag in all_frags.iter() {
        length_vec.push(frag.seq_dict.len() as u32);
    }
    length_vec.sort();
    return length_vec[(length_vec.len() as f64 * quantile) as usize];
}

pub fn get_length_gn(all_frags: &Vec<Frag>) -> SnpPosition {
    let mut last_pos = 0;
    let mut first_pos = std::u32::MAX;
    for frag in all_frags.iter() {
        if frag.last_position > last_pos {
            last_pos = frag.last_position;
        }
        if frag.first_position < first_pos {
            first_pos = frag.first_position;
        }
    }
    if last_pos > first_pos {
        return last_pos - first_pos;
    }
    else{
        return 0
    }
}

//Get the log p-value for a 1-sided binomial test. This is a asymptotically tight large deviation
//bound. It's super accurate when k/n >> p, but relatively inaccurate when k/n is close to p. One
//super nice thing about this approximation is that it is written as p = exp(A), so log(p) = A
//hence it is extremely numerically stable.
//
//I'm currently using this implementation. We can still mess around with using different approximations.
pub fn stable_binom_cdf_p_rev(n: usize, k: usize, p: f64, div_factor: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }

    //    return norm_approx(n,k,p,div_factor);

    let n64 = n as f64;
    let k64 = k as f64;

    //In this case, the relative entropy is bigger than the minimum of 0 which we don't want.
    let mut a = k64 / n64;

    if a == 1.0 {
        //Get a NaN error if a = 1.0;
        a = 0.9999999
    }
    if a == 0.0 {
        //Get a NaN error if we only have errors -- which can happen if we use polishing.
        a = 0.0000001;
    }
    

    let mut rel_ent = a * (a / p).ln() + (1.0 - a) * ((1.0 - a) / (1.0 - p)).ln();

    //If smaller error than epsilon, invert the rel-ent so that we get a positive probability
    //makes heuristic sense because a smaller than epsilon error is better than an epsilon error
    //for which the relative entropy is 0.
    if a < p {
        rel_ent = -rel_ent;
    }
    let large_dev_val = -1.0 * n64 / div_factor * rel_ent;
    //- 0.5 * (6.283*a*(1.0-a)*n64/div_factor).ln();

    return large_dev_val;

    //    return -1.0 * n64 / div_factor * rel_ent;
}

pub fn log_sum_exp(probs: &Vec<f64>) -> f64 {
    let max = probs.iter().copied().fold(f64::NAN, f64::max);
    let mut sum = 0.0;
    for logpval in probs {
        sum += (logpval - max).exp();
    }

    return max + sum.ln();
}

//Get a vector of read frequencies and error rates from a partition and its corresponding
//haplotype block.
pub fn get_seq_err_correlations(
    partition: &Vec<FxHashSet<&Frag>>,
    hap_block: &HapBlock,
    gap: SnpPosition,
) -> f64 {
    let ploidy = partition.len();
    let mut seq_err_corr_vec = vec![];
    for _i in 0..ploidy {
        let seq_err_corr_map_i = FxHashMap::default();
        seq_err_corr_vec.push(seq_err_corr_map_i);
    }
    for i in 0..ploidy {
        let haplo = &hap_block.blocks[i];
        for frag in partition[i].iter() {
            err_correlations(frag, haplo, &mut seq_err_corr_vec[i], gap);
        }
    }

    let mut list_of_p_values = vec![];
    let mut expected_diff_sum = 0.0;
    let mut max_diff_sum = 0.0;
    let mut _max_diff_sum_ploidy = 0.0;
    let mut max_diff_sum_chi = 0.0;
    let mut _pos_max_diff_sum_ploidy = 0;
    for j in 0..ploidy {
        for i in seq_err_corr_vec[j].keys() {
            let mut running_diff_sum = 0.0;
            let mut running_diff_sum_chi = 0.0;

            let mut errors_1 = FxHashMap::default();
            let mut errors_2 = FxHashMap::default();
            let mut total_alleles = 0.0;

            //            dbg!(&seq_err_corr_vec[j][i]);
            for (key, value) in seq_err_corr_vec[j].get(&i).unwrap().iter() {
                let e1 = errors_1.entry(key.0).or_insert(0.0);
                *e1 += *value as f64;
                total_alleles += *value as f64;
                let e2 = errors_2.entry(key.1).or_insert(0.0);
                *e2 += *value as f64;
            }
            for (_key, value) in errors_1.iter_mut() {
                *value /= total_alleles;
            }
            for (_key, value) in errors_2.iter_mut() {
                *value /= total_alleles;
            }

            let mut expected_map = FxHashMap::default();
            for (key, value) in errors_1.iter() {
                for (key2, value2) in errors_2.iter() {
                    expected_map.insert((*key, *key2), value * value2 * total_alleles);
                }
            }

            for (key, value) in expected_map {
                let diff_sum =
                    (value - (*seq_err_corr_vec[j][i].get(&key).unwrap_or(&0) as f64)).abs();
                let diff_sum_chi =
                    (value - (*seq_err_corr_vec[j][i].get(&key).unwrap_or(&0) as f64)).powf(2.0)
                        / value;

                running_diff_sum += diff_sum;
                running_diff_sum_chi += diff_sum_chi;
                expected_diff_sum += diff_sum;
            }
            if running_diff_sum > max_diff_sum {
                max_diff_sum = running_diff_sum;
            }
            if running_diff_sum_chi > max_diff_sum_chi {
                max_diff_sum_chi = running_diff_sum_chi;
            }
            let rv = ChiSquared::new(3.0).unwrap();
            let rv_res = 1.0 - rv.cdf(running_diff_sum_chi);
            list_of_p_values.push(rv_res);
        }
    }

    //    dbg!(
    //        expected_diff_sum,
    //        max_diff_sum_ploidy,
    //        pos_max_diff_sum_ploidy,
    //    );

    expected_diff_sum
}

pub fn err_correlations(
    r1: &Frag,
    hap: &Haplotype,
    seq_err_corr_map: &mut FxHashMap<SnpPosition, FxHashMap<(Genotype, Genotype), usize>>,
    gap: SnpPosition,
) {
    for pos in r1.positions.iter() {
        if !hap.contains_key(pos) {
            continue;
        }
        let current_var = *r1.seq_dict.get(pos).unwrap();

        if r1.seq_dict.contains_key(&(*pos + gap)) {
            let next_var = *r1.seq_dict.get(&(*pos + gap)).unwrap();
            //last_pos-1 because positions are 1-indexed
            let index = seq_err_corr_map
                .entry(*pos - 1)
                .or_insert(FxHashMap::default());
            let count = index.entry((current_var, next_var)).or_insert(0);
            *count += 1;
        }
    }
}

pub fn split_part_using_breaks<'a>(
    breaks: &FxHashMap<SnpPosition, FxHashSet<usize>>,
    part_to_split: &Vec<FxHashSet<&Frag>>,
    all_reads: &'a Vec<Frag>,
) -> Vec<Vec<FxHashSet<&'a Frag>>> {
    let mut breaks_with_min = FxHashMap::default();
    for (key, value) in breaks {
        if value.len() > 1 {
            breaks_with_min.insert(key, value);
        }
    }
    let ploidy = part_to_split.len();
    let mut split_parts = vec![vec![FxHashSet::default(); ploidy]; breaks_with_min.len() + 1];
    for (i, hap) in part_to_split.iter().enumerate() {
        for read in hap.iter() {
            if breaks_with_min.len() == 0 {
                split_parts[0][i].insert(&all_reads[read.counter_id]);
            } else {
                for (j, break_pos) in breaks_with_min.keys().sorted().enumerate() {
                    if read.first_position <= **break_pos {
                        split_parts[j][i].insert(&all_reads[read.counter_id]);
                    }
                    if read.last_position >= **break_pos {
                        split_parts[j + 1][i].insert(&all_reads[read.counter_id]);
                    }
                }
            }
        }
    }
    return split_parts;
}

pub fn get_range_with_lengths(
    snp_to_genome_pos: &Vec<usize>,
    block_length: usize,
    overlap_len: usize,
    minimal_density: f64,
) -> Vec<(SnpPosition, SnpPosition)> {
    let mut return_vec = vec![];
    let mut cum_pos = 0;
    let mut last_pos = snp_to_genome_pos[0];
    let mut left_endpoint = 0;
    let mut new_left_end = 0;
    let mut hit_new_left = false;

    for (i, pos) in snp_to_genome_pos.iter().enumerate() {
        let i = i as SnpPosition; //1-indexed
        if i == (snp_to_genome_pos.len() - 1) as SnpPosition {
            return_vec.push((left_endpoint, i));
            break;
        }
        if *pos < last_pos{
            log::error!("VCF malformed. Positions are not increasing {} {}", last_pos, *pos);
            std::process::exit(1);
        }
        cum_pos += *pos - last_pos;
        last_pos = *pos;
        if cum_pos > block_length - overlap_len && hit_new_left == false {
            new_left_end = i;
            hit_new_left = true;
        }
        if cum_pos > block_length {
            cum_pos = 0;
            let snp_density = (i - left_endpoint) as f64 / block_length as f64;
            if snp_density > minimal_density {
                return_vec.push((left_endpoint, i - 1));
            } else {
                log::debug!(
                    "Block endpoints {} - {} has density {} < min density {}",
                    left_endpoint,
                    i - 1,
                    snp_density,
                    minimal_density
                );
            }
            //left_endpoint = new_left_end;
            if snp_to_genome_pos[new_left_end as usize] + block_length
                < snp_to_genome_pos[(new_left_end + 1) as usize]
            {
                left_endpoint = new_left_end;
            } else {
                left_endpoint = new_left_end + 1;
            }
            last_pos = snp_to_genome_pos[left_endpoint as usize];
            hit_new_left = false;
        }
    }

    return_vec = return_vec.into_iter().map(|x| (x.0 + 1, x.1 + 1)).collect();
    return return_vec;
}

pub fn add_read_to_block(block: &mut HapBlock, frag: &Frag, part: usize) {
    for pos in frag.positions.iter() {
        let var_at_pos = frag.seq_dict.get(pos).unwrap();
        let sites = block.blocks[part]
            .entry(*pos)
            .or_insert(FxHashMap::default());
        let site_counter = sites.entry(*var_at_pos).or_insert(OrderedFloat(0.));
        *site_counter += phred_scale(frag, pos);
    }
}

pub fn remove_read_from_block(block: &mut HapBlock, frag: &Frag, part: usize) {
    for pos in frag.positions.iter() {
        let var_at_pos = frag.seq_dict.get(pos).unwrap();
        let sites = block.blocks[part]
            .entry(*pos)
            .or_insert(FxHashMap::default());
        let site_counter = sites.entry(*var_at_pos).or_insert(OrderedFloat(0.));
        if *site_counter != 0. {
            *site_counter -= phred_scale(frag, pos);
        }
        if *site_counter <= OrderedFloat(0.) {
            sites.remove(var_at_pos);
        }
    }
}

pub fn hybrid_correction(frags: Vec<Frag>) -> (Vec<Frag>, Vec<Frag>) {
    let start_t = Instant::now();
    let mut pos_to_frags = FxHashMap::default();
    let mut long_frags = vec![];

    for frag in frags.iter() {
        if frag.is_paired {
            for pos in frag.positions.iter() {
                let vec_pos = pos_to_frags.entry(pos).or_insert(FxHashSet::default());
                vec_pos.insert(frag);
            }
        } else {
            long_frags.push(frag);
        }
    }

    let final_frags: Mutex<Vec<_>> = Mutex::new(vec![]);
    long_frags.into_par_iter().for_each(|long_frag| {
        let mut covered_positions = FxHashSet::default();
        let mut covering_frags = FxHashSet::default();
        let mut positions = long_frag.positions.iter().collect::<Vec<&SnpPosition>>();
        positions.sort();
        for (i, _position) in positions.iter().enumerate() {
            if covered_positions.contains(positions[i]) {
                continue;
            }
            let mut covering_i_frags = FxHashSet::default();
            let mut j = i;
            loop {
                //Test
                if j != i {
                    break;
                }
                if j >= positions.len() {
                    break;
                }
                if !pos_to_frags.contains_key(&positions[j]) {
                    break;
                }
                let covering_i = &pos_to_frags[&positions[j]];
                if covering_i.is_disjoint(&covering_i_frags) && j != i {
                    break;
                }
                if j == i {
                    covering_i_frags = covering_i_frags.union(&covering_i).copied().collect();
                } else {
                    covering_i_frags = covering_i_frags
                        .intersection(&covering_i)
                        .copied()
                        .collect();
                }
                j += 1;
            }
            if covering_i_frags.is_empty() {
                continue;
            }
            let best_frag = covering_i_frags
                .into_iter()
                .max_by_key(|x| {
                    let d = distance(x, &long_frag);
                    (d.0 * 10) / (d.1 + 1)
                })
                .unwrap();
            for pos in best_frag.positions.iter() {
                covered_positions.insert(*pos);
            }
            covering_frags.insert(best_frag);
        }
        let cand_seq_dict = set_to_seq_dict(&covering_frags, true);
        let mut locked = final_frags.lock().unwrap();
        locked.push(correct_long_read(&cand_seq_dict, &long_frag));
        //        dbg!(cand_seq_dict, &long_frag.seq_dict);
    });

    log::info!("Time taken error_correct {:?}", Instant::now() - start_t);
    let mut short_frags = vec![];
    for frag in frags.into_iter() {
        if frag.is_paired {
            short_frags.push(frag);
        }
    }
    return (final_frags.into_inner().unwrap(), short_frags);
}

fn correct_long_read(short_frags_dict: &Haplotype, long_frag: &Frag) -> Frag {
    let mut new_frag = long_frag.clone();
    for pos in new_frag.positions.iter() {
        if !short_frags_dict.contains_key(pos) {
            continue;
        }
        let val = new_frag.seq_dict.get_mut(pos).unwrap();
        if short_frags_dict[&pos].len() > 1 {
            continue;
        } else {
            *val = *short_frags_dict[&pos]
                .iter()
                .max_by_key(|entry| entry.1)
                .unwrap()
                .0;
        }
    }
    return new_frag;
}

pub fn get_errors_cov_from_frags(
    frags: &FxHashSet<&Frag>,
    left_snp_pos: SnpPosition,
    right_snp_pos: SnpPosition,
) -> (f64, f64, f64, f64) {
    let mean = true;
    //don't use phred here while computing the SNP error rate
    //for outputs.
    let hap_map = set_to_seq_dict(frags, false);
    let mut snp_counter_list = vec![];
    let emptydict = FxHashMap::default();
    let mut errors = 0.;
    let mut total_support = 0.;
    let mut snp_nonzero = FxHashSet::default();
    for pos in left_snp_pos..right_snp_pos + 1 {
        let mut snp_support = OrderedFloat(0.);
        let mut max_count_pos = OrderedFloat(0.);
        let allele_map = hap_map.get(&pos).unwrap_or(&emptydict);
        if *allele_map != emptydict {
            snp_nonzero.insert(pos);
            for (site, count) in allele_map {
                if *site == GAP_CHAR {
                    continue;
                }
                if *count > snp_support {
                    max_count_pos = *count;
                }
                snp_support += count;
            }
        }
        total_support += snp_support.into_inner();
        errors += snp_support.into_inner() - max_count_pos.into_inner();
        snp_counter_list.push(snp_support);
    }
    snp_counter_list.sort();
    let cov;
    if snp_counter_list.is_empty() {
        cov = 0.;
    } else {
        //Quantile
        if !mean {
            cov = snp_counter_list[snp_counter_list.len() * 2 / 3].into_inner();
        }
        //Mean
        else {
            //cov = *snp_counter_list.iter().sum::<GenotypeCount>() / snp_counter_list.len() as f64;
            if snp_nonzero.len() > 0{
                cov = *snp_counter_list.iter().sum::<GenotypeCount>() / snp_nonzero.len() as f64;
            }
            else{
                cov = 0.;
            }
        }
    }

    return (
        cov,
        errors as f64 / total_support as f64,
        errors as f64,
        total_support as f64,
    );
}

fn get_iqr(hap: &Haplotype) -> [OrderedFloat<f64>; 2] {
    let ret;
    let mut total_covs = vec![];
    for (_pos, geno_dict) in hap.iter() {
        let total_cov = geno_dict.iter().map(|x| *x.1).sum::<GenotypeCount>();
        total_covs.push(total_cov);
    }
    if total_covs.is_empty() {
        ret = [OrderedFloat(0.), OrderedFloat(f64::MAX/2.)];
    } else {
        total_covs.sort();
        let q1 = total_covs[total_covs.len() / 4];
        let q4 = total_covs[total_covs.len() * 3 / 4];
        let iqr = q4 - q1;
        ret = [q1 - iqr * OrderedFloat(2.), q4 + OrderedFloat(2.) * iqr];
    }
    return ret;
}

pub fn distance_between_haplotypes(
    hap1: &Haplotype,
    hap2: &Haplotype,
    range: &(SnpPosition, SnpPosition),
    reliability_cutoff: f64,
    cov_cutoff: f64,
) -> (f64, f64) {
    let hap1_iqr_intervals = get_iqr(hap1);
    let hap2_iqr_intervals = get_iqr(hap2);
    log::trace!("hap1_iqr_intervals {:?}", hap1_iqr_intervals);
    log::trace!("hap2_iqr_intervals {:?}", hap2_iqr_intervals);

    let mut diff_poses = vec![];
    let cov_cutoff = OrderedFloat(cov_cutoff);
    let mut same = 0.;
    let mut diff = 0.;
    let ratioshap1 = hap1.iter().map(|(snp, genodict)| {
        let total_cov = genodict.iter().map(|x| *x.1).sum::<GenotypeCount>();
        let highest_cov = genodict.iter().map(|x| *x.1).max().unwrap();
        let ratio = highest_cov / total_cov;
        (snp, ratio)
    }).collect::<FxHashMap<_, _>>();
    let ratioshap2 = hap2.iter().map(|(snp, genodict)| {
        let total_cov = genodict.iter().map(|x| *x.1).sum::<GenotypeCount>();
        let highest_cov = genodict.iter().map(|x| *x.1).max().unwrap();
        let ratio = highest_cov / total_cov;
        (snp, ratio)
    }).collect::<FxHashMap<_, _>>();
    for pos in hap1.keys() {
        let cov_pos_1 = hap1[pos].iter().map(|x| *x.1).sum::<GenotypeCount>();
        if hap2.contains_key(pos) {

            if ratioshap1[pos].into_inner() < reliability_cutoff || ratioshap2[pos].into_inner() < reliability_cutoff {
                continue;
            }
            let cov_pos_2 = hap2[pos].iter().map(|x| *x.1).sum::<GenotypeCount>();

            let pos_hap1_fits_iqr = cov_pos_1 > hap1_iqr_intervals[0]
                && cov_pos_1 < hap1_iqr_intervals[1];
            let pos_hap2_fits_iqr = cov_pos_2 > hap2_iqr_intervals[0]
                && cov_pos_2 < hap2_iqr_intervals[1];

            if (cov_pos_1 > cov_cutoff && cov_pos_2 > cov_cutoff)
                && (*pos >= range.0 && *pos <= range.1)
                && pos_hap1_fits_iqr
                && pos_hap2_fits_iqr
            {
                let consensus_var1 = hap1
                    .get(pos)
                    .unwrap()
                    .iter()
                    .max_by_key(|entry| entry.1)
                    .unwrap()
                    .0;

                let consensus_var2 = hap2
                    .get(pos)
                    .unwrap()
                    .iter()
                    .max_by_key(|entry| entry.1)
                    .unwrap()
                    .0;

                if consensus_var1 == consensus_var2 {
                    //same += (((ratioshap1[pos] - 0.5) * (ratioshap2[pos] - 0.5))/(0.25)).into_inner();
                    same += 1.
                } else {
                    //diff += (((ratioshap1[pos] - 0.5) * (ratioshap2[pos] - 0.5))/(0.25)).into_inner();
                    diff += 1.;
                    diff_poses.push(*pos);
                }
            }
        }
    }

    log::trace!("diff_poses {:?}", diff_poses);

    return (same, diff);
}

pub fn phred_scale_dbg(frag: &FragDBG, pos: &SnpPosition) -> GenotypeCount {
    if constants::USE_QUAL_SCORES {
        let qual_score = frag.qual_dict.get(&pos).unwrap();
        let prob = 1. - 10_f32.powf((*qual_score) as f32 / -10.);
        return OrderedFloat(prob.into());
    } else {
        return OrderedFloat(1.);
    }
}

#[inline]
pub fn phred_scale(frag: &Frag, pos: &SnpPosition) -> GenotypeCount {
    if constants::USE_QUAL_SCORES {
        let qual_score = frag.qual_dict.get(&pos).unwrap();
        let prob = 1. - 10_f32.powf((*qual_score) as f32 / -10.);
        return OrderedFloat(prob.into());
    } else {
        return OrderedFloat(1.);
    }
}

pub fn remove_monomorphic_allele(mut frags: Vec<Frag>, error: f64) -> Vec<Frag> {
    let mut allele_count_map: FxHashMap<SnpPosition, FxHashMap<Genotype, GenotypeCount>> =
        FxHashMap::default();
    let mut mono_alleles: FxHashSet<SnpPosition> = FxHashSet::default();
    let mut new_frags = vec![];

    for frag in frags.iter() {
        for (snp_pos, geno) in frag.seq_dict.iter() {
            let count_m = allele_count_map
                .entry(*snp_pos)
                .or_insert(FxHashMap::default());
            let count = count_m.entry(*geno).or_insert(OrderedFloat(0.));
            *count += phred_scale(frag, snp_pos);
        }
    }

    for (allele, map) in allele_count_map.iter() {
        if map.len() == 1 {
            mono_alleles.insert(*allele);
            log::trace!("allele {} removed", allele);
        } else {
            let mut vals = map.values().collect::<Vec<&GenotypeCount>>();
            vals.sort_by(|x, y| y.partial_cmp(&x).unwrap());
            if vals[0].into_inner() * error > vals[1].into_inner() {
                mono_alleles.insert(*allele);
                log::trace!("allele {} removed, {:?}", allele, map);
            }
        }
    }

    for frag in frags.iter_mut() {
        let mut to_rem = vec![];
        for pos in frag.positions.iter() {
            if mono_alleles.contains(pos) {
                frag.seq_dict.remove(pos);
                frag.qual_dict.remove(pos);
                frag.snp_pos_to_seq_pos.remove(pos);
                to_rem.push(*pos);
            }
        }
        for pos in to_rem {
            frag.positions.remove(&pos);
        }
        let mut positions_new = frag.positions.iter().collect::<Vec<&SnpPosition>>();
        positions_new.sort();
        if !positions_new.is_empty() {
            frag.first_position = **positions_new.first().unwrap();
            frag.last_position = **positions_new.last().unwrap();
            new_frags.push(mem::take(frag));
        }
    }

    new_frags.sort();
    for (i, frag) in new_frags.iter_mut().enumerate() {
        frag.counter_id = i;
    }

    return new_frags;

}
