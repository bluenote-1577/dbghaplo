use crate::types_structs::*;
use fxhash::{FxHashMap, FxHashSet};
use petgraph::dot::{Config, Dot};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::{algo, prelude::*};
use statrs::distribution::{Binomial, Discrete, DiscreteCDF};
use statrs::statistics::Distribution;
use std::collections::VecDeque;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::io::{self, BufReader, BufWriter, Write};
use std::mem;
use std::rc::Rc;

pub fn frag_to_dbgfrag(frag: &Frag) -> FragDBG {
    let mut seq = vec![];
    for (pos, geno) in frag.seq_dict.iter() {
        seq.push((*pos, *geno));
    }
    seq.sort_by(|a, b| a.0.cmp(&b.0));
    let toret = FragDBG {
        id: frag.id.clone(),
        counter_id: frag.counter_id,
        seq,
        seq_dict: frag.seq_dict.clone(),
        seq_string: frag.seq_string.clone(),
        first_position: frag.first_position,
        last_position: frag.last_position,
        snp_pos_to_seq_pos: frag.snp_pos_to_seq_pos.clone(),
        first_pos_base: frag.first_pos_base,
        last_pos_base: frag.last_pos_base,
    };
    toret
}

#[derive(Debug, Clone, Eq, PartialEq, Default)]
pub struct DBGInfo {
    pub out_varmers: Vec<Rc<VarMer>>,
    pub in_varmers: Vec<Rc<VarMer>>,
    pub coverage: u64,
}

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
struct VarmerPath {
    pub varmers: Vec<DictFrag>,
    pub total_avg_cov: u64,
}

struct VarmerPathInteger {
    pub intver: VarMer,
    pub total_avg_cov: u64,
}

impl fmt::Display for VarmerPath {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let seq1 = self
            .varmers
            .iter()
            .map(|varmer| fmt_varmer_d(varmer))
            .collect::<Vec<String>>()
            .join("->");
        write!(f, "{}", seq1)
    }
}

fn get_edges_varmers(dbg: &mut FxHashMap<VarMer, DBGInfo>, k: usize) {
    let mut suffixes: FxHashMap<VarMer, Vec<Rc<VarMer>>> = FxHashMap::default();
    let mut prefixes: FxHashMap<VarMer, Vec<Rc<VarMer>>> = FxHashMap::default();

    for (varmer, _info) in dbg.iter() {
        //dbg!(varmer.len(),k);
        let suffix = varmer[varmer.len() - k + 1..].to_vec();
        let prefix = varmer[..k - 1].to_vec();
        let varmer_ref = Rc::new(varmer.clone());

        suffixes
            .entry(suffix)
            .or_insert_with(Vec::new)
            .push(varmer_ref.clone());
        prefixes
            .entry(prefix)
            .or_insert_with(Vec::new)
            .push(varmer_ref.clone());
    }

    for (prefix, out_varmers) in prefixes.iter() {
        if let Some(in_varmers) = suffixes.get(prefix) {
            for in_varmer in in_varmers.iter() {
                for out_varmer in out_varmers.iter() {
                    if let Some(in_info) = dbg.get_mut(&**in_varmer) {
                        in_info.out_varmers.push(Rc::clone(out_varmer));
                    }
                    if let Some(out_info) = dbg.get_mut(&**out_varmer) {
                        out_info.in_varmers.push(Rc::clone(in_varmer));
                    }
                }
            }
        }
    }
}

pub fn print_dbg(dbg: &FxHashMap<VarMer, DBGInfo>, file_name: &str) {
    let mut dot = String::from("digraph G {\n");
    for (node, info) in dbg.iter() {
        let out_edges = &info.out_varmers;
        if out_edges.len() == 0 && info.in_varmers.len() == 0 {
            let mut s1 = node
                .iter()
                .map(|(_, geno)| format!("{}", geno))
                .collect::<Vec<String>>()
                .join("");
            s1.push_str(&format!(
                "_{}-{}:{}-LEN:{}",
                node.first().unwrap().0,
                node.last().unwrap().0,
                info.coverage,
                node.len()
            ));
            dot.push_str(&format!("    \"{}\" -> \"{}\";\n", s1, s1));
        }
        for out in out_edges {
            let mut s1 = node
                .iter()
                .map(|(_, geno)| format!("{}", geno))
                .collect::<Vec<String>>()
                .join("");
            let mut s2 = out
                .iter()
                .map(|(_, geno)| format!("{}", geno))
                .collect::<Vec<String>>()
                .join("");
            //append coverage to s1 and s2
            s1.push_str(&format!(
                "_{}-{}:{}-LEN:{}",
                node.first().unwrap().0,
                node.last().unwrap().0,
                info.coverage,
                node.len()
            ));
            s2.push_str(&format!(
                "_{}-{}:{}-LEN:{}",
                out.first().unwrap().0,
                out.last().unwrap().0,
                dbg.get(&**out).unwrap().coverage,
                out.len()
            ));
            dot.push_str(&format!("    \"{}\" -> \"{}\";\n", s1, s2));
        }
    }
    dot.push_str("}\n");

    //write to file called t.dot
    std::fs::write(file_name, dot).expect("Unable to write file");
}

fn dbg_from_frags(
    dbg_frags: &Vec<FragDBG>,
    k: usize,
    prev_dbg: Option<FxHashMap<VarMer, DBGInfo>>,
    unitigs: Option<&FxHashMap<VarMer, DBGInfo>>,
    base_k: Option<usize>,
) -> FxHashMap<VarMer, DBGInfo> {
    // DBG Format: key -> edges
    let mut dbg: FxHashMap<VarMer, DBGInfo> = FxHashMap::default();
    let mut sequences = vec![];

    if let Some(unitigs) = unitigs.as_ref() {
        for (unitig, _) in unitigs.iter() {
            sequences.push(unitig);
        }
    }
    for frag in dbg_frags.iter() {
        sequences.push(&frag.seq);
    }

    for seq in sequences.iter() {
        let mut varmer = VecDeque::new();
        for (curr_pos, curr_geno) in seq.into_iter() {
            if *curr_pos > 250 {
                break;
            }
            varmer.push_back((*curr_pos, *curr_geno));
            if varmer.len() == k {
                let new_varmer = varmer.iter().cloned().collect::<VarMer>();
                if let Some(prev_dbg) = prev_dbg.as_ref() {
                    let base_k = base_k.unwrap();
                    let prev_kmer = &new_varmer[0..k - base_k];
                    let end_base_kmer = &new_varmer[base_k..];
                    let min_cov;
                    if let Some(prev_info) = prev_dbg.get(prev_kmer) {
                        if let Some(end_cov) = prev_dbg.get(end_base_kmer) {
                            min_cov = prev_info.coverage.min(end_cov.coverage);
                        } else {
                            min_cov = 0;
                        }
                    } else {
                        min_cov = 0;
                    }

                    if min_cov != 0 {
                        dbg.entry(new_varmer).or_insert(DBGInfo {
                            out_varmers: vec![],
                            in_varmers: vec![],
                            coverage: min_cov,
                        });
                    }
                } else {
                    let info = dbg.entry(new_varmer).or_insert(DBGInfo {
                        out_varmers: vec![],
                        in_varmers: vec![],
                        coverage: 0,
                    });
                    info.coverage += 1;
                }
                varmer.pop_front();
            }
        }
    }

    get_edges_varmers(&mut dbg, k);
    return dbg;
}

fn filter_dbg(
    dbg: FxHashMap<VarMer, DBGInfo>,
    min_cov: Option<u64>,
    bad_nodes: Option<Vec<VarMer>>,
    k: usize,
) -> FxHashMap<VarMer, DBGInfo> {
    let mut new_dbg = FxHashMap::default();
    for (node, info) in dbg.iter() {
        if let Some(bad_nodes) = bad_nodes.as_ref() {
            if bad_nodes.contains(node) {
                continue;
            }
        }
        if let Some(min_cov) = min_cov {
            if info.coverage <= min_cov {
                continue;
            }
        }
        new_dbg.insert(
            node.clone(),
            DBGInfo {
                out_varmers: vec![],
                in_varmers: vec![],
                coverage: info.coverage,
            },
        );
    }
    get_edges_varmers(&mut new_dbg, k);
    return new_dbg;
}

fn get_unitigs(
    dbg: &FxHashMap<VarMer, DBGInfo>,
    k: usize,
    unitig_graph: bool,
) -> FxHashMap<VarMer, DBGInfo> {
    let mut unitigs = FxHashMap::default();
    let mut visited_nodes = FxHashSet::default();
    for (node, info) in dbg.iter() {
        if visited_nodes.contains(node) {
            continue;
        }
        if info.in_varmers.len() != 1
            || (info.in_varmers.len() == 1
                && dbg.get(&**info.in_varmers[0]).unwrap().out_varmers.len() != 1)
        {
            let mut cov = info.coverage;
            let mut count = 1;
            visited_nodes.insert(node);
            let mut unitig = node.clone();

            if info.out_varmers.len() == 1 {
                let mut prev = &node.clone();
                let mut next = &info.out_varmers[0];
                while dbg.get(&**prev).unwrap().out_varmers.len() == 1
                    && dbg.get(&**next).unwrap().in_varmers.len() == 1
                {
                    visited_nodes.insert(next);
                    if unitig_graph {
                        for pg in next.iter() {
                            if unitig.contains(pg) {
                                continue;
                            } else {
                                unitig.push(*pg);
                            }
                        }
                    } else {
                        unitig.push(next.last().unwrap().clone());
                    }

                    cov += dbg.get(&**next).unwrap().coverage;
                    count += 1;

                    if dbg.get(&**next).unwrap().out_varmers.len() == 0 {
                        break;
                    }
                    prev = next;
                    next = &dbg.get(&**next).unwrap().out_varmers[0];
                }
            }
            unitigs.insert(
                unitig,
                DBGInfo {
                    out_varmers: vec![],
                    in_varmers: vec![],
                    coverage: cov / count,
                },
            );
        }
    }
    get_edges_varmers(&mut unitigs, k);
    return unitigs;
}

fn count_kmers(dbg_frags: &Vec<FragDBG>, k: usize) -> (FxHashMap<VarMer, u64>, FxHashSet<u32>) {
    let mut snps = FxHashSet::default();
    let mut kmers = FxHashMap::default();
    for frag in dbg_frags.iter() {
        let seq = &frag.seq;
        let mut varmer = VecDeque::new();
        for (curr_pos, curr_geno) in seq.iter() {
            snps.insert(*curr_pos);
            if *curr_pos > 250 {
                break;
            }
            varmer.push_back((*curr_pos, *curr_geno));
            if varmer.len() == k {
                let new_varmer = varmer.iter().cloned().collect::<VarMer>();
                *kmers.entry(new_varmer).or_insert(0) += 1;
                varmer.pop_front();
            }
        }
    }
    return (kmers, snps);
}

pub fn construct_dbg(dbg_frags: &Vec<FragDBG>, k: usize, end: usize) {
    let (_, snps) = count_kmers(dbg_frags, k);
    let k = k.min(snps.len() * 7 / 10).max(1);
    let end = (end).min(snps.len() * 8 / 10).max(k);
    let end = end - k;
    log::info!("Start k: {}, end k: {}", k, end);

    let dbg = dbg_from_frags(dbg_frags, k, None, None, None);
    let (mut kmer_count, snps) = count_kmers(dbg_frags, k);
    let total_cov = kmer_count.iter().fold(0, |acc, (_varmer, cov)| acc + cov);
    let min_cov = u64::max(total_cov / (snps.len() as u64 - k as u64 + 1) / 200, 2);
    log::info!("Min Cov: {:?}", min_cov);

    let mut dbg = filter_dbg(dbg, Some(min_cov), None, k);
    let mut uni = get_unitigs(&dbg, k, false);
    kmer_count.retain(|varmer, _cov| dbg.contains_key(varmer));
    let step = 1;

    for l in (step..end + 1).step_by(step) {
        dbg = dbg_from_frags(dbg_frags, l + k, Some(dbg), Some(&uni), Some(step));
        uni = get_unitigs(&dbg, k + l, false);
    }

    print_dbg(&dbg, "dbg.dot");
    //Unitigging
    let unitigs = uni;
    print_dbg(&unitigs, "unitigs.dot");

    let mut final_unitigs = unitigs;
    for i in 1..3 {
        let bad_unitigs = query_unitigs(&final_unitigs, i);
        log::info!("Number of bad unitigs {}", bad_unitigs.len());
        let filtered_unitigs = filter_dbg(final_unitigs, None, Some(bad_unitigs), k + end);
        final_unitigs = get_unitigs(&filtered_unitigs, k + end, true);
        print_dbg(&final_unitigs, "clean_unitigs.dot");
    }

    let final_unitigs = clean_hanging_kmers(final_unitigs, k + end - 1);
    print_dbg(&final_unitigs, "nohang_unitigs.dot");

    //Try aligning reads to graph
    let vec_df_all: Vec<DictFrag> = dbg_frags
        .iter()
        .map(|frag| fragdbg_to_dictfrag(frag))
        .collect();
    let vec_df_graph: Vec<DictFrag> = final_unitigs
        .iter()
        .map(|(varmer, info)| DictFrag {
            seq: varmer.iter().cloned().collect(),
            seq_vec: varmer.clone(),
            first_position: varmer.first().unwrap().0,
            last_position: varmer.last().unwrap().0,
            cov: info.coverage,
        })
        .collect();

    let dict_frag_to_index = vec_df_graph
        .iter()
        .enumerate()
        .map(|(i, frag)| (frag, i))
        .collect::<FxHashMap<&DictFrag, usize>>();

    let mut path_dict = FxHashMap::default();
    for dict_frag in vec_df_all.iter() {
        log::info!("ALIGN TO GRAPH");
        print_varmer_d(dict_frag);
        let hits = get_hits(dict_frag, &vec_df_graph, 100, true);
        let dp_res = dp_hits(&hits, dict_frag, &final_unitigs, 100, -1.0, false);
        let varmers = varmers_from_dp_res(&dp_res);
        log::info!("HITS: {:?}", hits.len());
        log::info!("DP RES: {:?}", dp_res.score);
        for varmer in varmers.iter() {
            print_varmer_d(varmer);
        }
        *path_dict.entry(varmers).or_insert(0) += 1;
    }

    let mut unitig_paths = vec![];

    for (varmers, count) in path_dict {
        log::info!("PATH COUNT: {}", count);

        for varmer in varmers.iter() {
            print_varmer_d(varmer);
        }
        let path = VarmerPath {
            varmers,
            total_avg_cov: count,
        };
        unitig_paths.push(path);
    }

    //1) collapse contained paths with coverage

    let mut graph_dict = FxHashMap::default();
    for path in unitig_paths.iter() {
        graph_dict.entry(path).or_insert((vec![], vec![]));
    }
    for i in 0..unitig_paths.len() {
        for j in 0..unitig_paths.len() {
            if i == j {
                continue;
            }
            let path1 = &unitig_paths[i];
            let path2 = &unitig_paths[j];
            let mut contained = false;

            if path1
                .varmers
                .iter()
                .collect::<FxHashSet<_>>()
                .is_subset(&path2.varmers.iter().collect::<FxHashSet<_>>())
            {
                contained = true;
            }

            if contained {
                let (out_edges, _) = graph_dict.get_mut(path1).unwrap();
                out_edges.push(path2);
                let (_, in_edges) = graph_dict.get_mut(path2).unwrap();
                in_edges.push(path1);
            }
        }
    }

    let mut outside_paths = FxHashMap::default();
    for path in unitig_paths.iter() {
        let (out_edges, _) = graph_dict.get(path).unwrap();
        if out_edges.len() == 0 {
            outside_paths.insert(path.varmers.clone(), path.total_avg_cov);
        }
    }

   for path in unitig_paths.iter(){
       let (out_edges, _) = graph_dict.get(&path).unwrap();
       let mut outside_out_edges = vec![];
       for path in out_edges.iter(){
           if outside_paths.contains_key(&path.varmers){
               outside_out_edges.push(path);
           }
       }
       if outside_out_edges.len() != 1{
           continue;
       }
       for outside_path in outside_out_edges.iter(){
           let len = outside_out_edges.len() as u64;
           let cov = outside_paths.get_mut(&outside_path.varmers).unwrap();
           *cov += path.total_avg_cov / len;
       }
   }

    for path in outside_paths.iter() {
        println!("OUTSIDE PATH: COV {}", outside_paths[path.0]);
        for varmer in path.0.iter() {
            print_varmer_d(varmer);
        }
    }

    let outside_paths = outside_paths
        .into_iter()
        .map(|(varmers, cov)| VarmerPath {
            varmers,
            total_avg_cov: cov,
        })
        .collect::<Vec<VarmerPath>>();
    let mut integer_paths = vec![];
    for path in outside_paths.iter() {
        let integer_node = path
            .varmers
            .iter()
            .map(|varmer| (0, *dict_frag_to_index.get(varmer).unwrap() as u8))
            .collect::<Vec<(u32, u8)>>();
        let integer_path = VarmerPathInteger {
            intver: integer_node,
            total_avg_cov: path.total_avg_cov,
        };
        integer_paths.push(integer_path);
    }

    //3) create overlap path graph

    let mut assembly_graph: FxHashMap<VarMer, DBGInfo> = FxHashMap::default();
    for path1 in integer_paths.iter() {
        assembly_graph
            .entry(path1.intver.clone())
            .or_insert(DBGInfo {
                out_varmers: vec![],
                in_varmers: vec![],
                coverage: path1.total_avg_cov,
            });
    }
    for path1 in integer_paths.iter() {
        for path2 in integer_paths.iter() {
            //check for suffix prefix overlaps of VarmerPaths.varmers
            let mut overlap_len = 0;
            for k in 0..path1.intver.len().min(path2.intver.len()) {
                if path1.intver[k..] == path2.intver[..path2.intver.len() - k] {
                    overlap_len = k;
                    break;
                }
            }
            if overlap_len > 1 {
                let info1 = assembly_graph.get_mut(&path1.intver).unwrap();
                info1.out_varmers.push(Rc::new(path2.intver.clone()));
                let info2 = assembly_graph.get_mut(&path2.intver).unwrap();
                info2.in_varmers.push(Rc::new(path1.intver.clone()));
            }
        }
    }

    print_dbg(&assembly_graph, "assembly_graph.dot");

    let integer_unitigs = get_unitigs(&assembly_graph, 1, true);
    let mut paths = vec![];
    for (int_unitig, info) in integer_unitigs.iter(){
        let path_as_df = int_unitig
            .iter()
            .map(|(_, x)| &vec_df_graph[(*x) as usize])
            .collect::<Vec<_>>();
        let min_cov_path = path_as_df.iter().map(|x| x.cov).min().unwrap();
        if min_cov_path > info.coverage * 2 {
            log::debug!(
                "INTEGER UNITIG FAILED -- MIN COV PATH: {} INFO COV: {}",
                min_cov_path,
                info.coverage
            );
            for df in path_as_df.iter() {
                print_varmer_d(df);
            }
        } else {
            log::debug!(
                "INTEGER UNITIG PASSED -- MIN COV PATH: {} INFO COV: {}",
                min_cov_path,
                info.coverage
            );
            for df in path_as_df.iter() {
                print_varmer_d(df);
            }
            let path_as_varmers = int_unitig
                .iter()
                .rev()
                .map(|(_, x)| (&vec_df_graph[(*x) as usize].seq_vec, info.coverage as usize))
                .collect::<Vec<_>>();
            paths.push(path_as_varmers);

        }
    }
    

    get_path_haps(dbg_frags, &final_unitigs, paths);
}

pub fn get_paths(
    unitigs: &FxHashMap<VarMer, DBGInfo>,
    conservative: bool,
) -> Vec<Vec<(&VarMer, usize)>> {
    let mut visited_nodes = FxHashSet::default();
    let mut paths = vec![];
    let mut sinks_and_forks = vec![];
    for (node, info) in unitigs.iter() {
        if info.in_varmers.len() == 0 {
            sinks_and_forks.push((node, info));
        } else if info.out_varmers.len() > 1 {
            sinks_and_forks.push((node, info));
        }
    }
    sinks_and_forks.sort_by(|a, b| b.1.coverage.cmp(&a.1.coverage));
    for (node, info) in sinks_and_forks.iter() {
        if info.out_varmers.len() == 0 {
            paths.push(vec![(*node, info.coverage as usize)]);
        }
        for out in info.out_varmers.iter() {
            let mut path = vec![(*node, info.coverage as usize)];
            let mut curr_node = &**out;
            loop {
                if visited_nodes.contains(curr_node) && conservative {
                    break;
                }
                visited_nodes.insert(curr_node.clone());
                let info = unitigs.get(&**curr_node).unwrap();
                path.push((curr_node, info.coverage as usize));
                let num_in_edges = info.in_varmers.len();
                if num_in_edges == 1 {
                    curr_node = &unitigs.get(&**curr_node).unwrap().out_varmers[0];
                } else {
                    paths.push(path);
                    break;
                }
            }
        }
    }
    return paths;
}

fn get_hits<'a>(
    varmer_d: &'a DictFrag,
    dict_frags: &'a Vec<DictFrag>,
    hit_threshold: usize,
    allow_equal: bool,
) -> Vec<Hit<'a>> {
    let mut hits = vec![];
    for ref_varmer in dict_frags.iter() {
        if varmer_d == ref_varmer && !allow_equal {
            continue;
        }

        let mut same = FxHashSet::default();
        let mut r_to_a = FxHashSet::default();
        let mut a_to_r = FxHashSet::default();
        let mut del = FxHashSet::default();

        //check if unitigs overlap based on first and last position
        if (varmer_d.first_position <= ref_varmer.last_position
            && varmer_d.last_position >= ref_varmer.first_position)
            || (ref_varmer.first_position <= varmer_d.last_position
                && ref_varmer.last_position >= varmer_d.first_position)
        {
            for (pos, geno) in ref_varmer.seq.iter() {
                if varmer_d.seq.contains_key(pos) {
                    let gn2 = *geno;
                    let gn1 = *varmer_d.seq.get(pos).unwrap();
                    if gn1 == gn2 {
                        same.insert(*pos);
                    } else if gn1 == 0 && gn2 != 0 {
                        r_to_a.insert(*pos);
                    } else if gn1 != 0 && gn2 == 0 {
                        a_to_r.insert(*pos);
                    } else {
                        a_to_r.insert(*pos);
                    }
                } else if *pos <= varmer_d.last_position && *pos >= varmer_d.first_position {
                    del.insert(*pos);
                }
            }

            for (pos, _) in varmer_d.seq.iter() {
                if !ref_varmer.seq.contains_key(pos)
                    && *pos <= ref_varmer.last_position
                    && *pos >= ref_varmer.first_position
                {
                    del.insert(*pos);
                }
            }
            if r_to_a.len() + a_to_r.len() + del.len() <= hit_threshold {
                hits.push(Hit {
                    varmer: &ref_varmer,
                    same,
                    r_to_a,
                    a_to_r,
                    del,
                });
            }
        }
    }

    hits.sort_by(|a, b| {
        (a.varmer.first_position, a.varmer.last_position)
            .cmp(&(b.varmer.first_position, b.varmer.last_position))
    });
    return hits;
}

fn varmers_from_dp_res<'a>(dp_res: &'a DpResult) -> Vec<DictFrag> {
    if dp_res.traceback.len() == 0 {
        return vec![];
    }
    let mut varmers = vec![];
    let mut traceback = dp_res.max_index;
    while traceback != dp_res.traceback[traceback] {
        varmers.push(dp_res.hits[traceback].varmer.clone());
        traceback = dp_res.traceback[traceback];
    }
    varmers.push(dp_res.hits[traceback].varmer.clone());
    return varmers;
}

fn dp_hits<'a>(
    hits: &'a Vec<Hit>,
    varmer_d: &DictFrag,
    graph: &FxHashMap<VarMer, DBGInfo>,
    threshold: usize,
    mismatch_pen: f64,
    ambiguous_allowed: bool,
) -> DpResult<'a> {


    let empty_result =  DpResult {
            score: (0., 0),
            traceback: vec![],
            max_index: 0,
            hits: &hits,
            dels_max: FxHashSet::default(),
            rtoa_max: FxHashSet::default(),
            ator_max: FxHashSet::default(),
            same_max: FxHashSet::default(),
            total_errs: 0,
        };

    if hits.len() == 0 {
        return empty_result
    }

    let mut dp_vec = hits
        .iter()
        .map(|x| {
            (
                x.same.len() as f64 + mismatch_pen * (x.r_to_a.len() + x.a_to_r.len() + x.del.len()) as f64,
                x.varmer.cov,
            )
        })
        .collect::<Vec<(f64, u64)>>();
    let mut traceback_vec = (0..hits.len()).map(|x| x).collect::<Vec<usize>>();
    //        let mut same_vecs : Vec<FxHashSet<u32>> = vec![FxHashSet::default(); hits.len()];
    //        let mut rtoa_vecs : Vec<FxHashSet<u32>> = vec![FxHashSet::default(); hits.len()];
    //        let mut ator_vecs : Vec<FxHashSet<u32>> = vec![FxHashSet::default(); hits.len()];
    //        let mut del_vecs : Vec<FxHashSet<u32>> = vec![FxHashSet::default(); hits.len()];
    let mut same_vecs = hits
        .iter()
        .map(|x| x.same.clone())
        .collect::<Vec<FxHashSet<u32>>>();
    let mut rtoa_vecs = hits
        .iter()
        .map(|x| x.r_to_a.clone())
        .collect::<Vec<FxHashSet<u32>>>();
    let mut ator_vecs = hits
        .iter()
        .map(|x| x.a_to_r.clone())
        .collect::<Vec<FxHashSet<u32>>>();
    let mut del_vecs = hits
        .iter()
        .map(|x| x.del.clone())
        .collect::<Vec<FxHashSet<u32>>>();

    for i in 0..hits.len() {
        let mut best_index = i;
        //let mut new_scores = vec![0.; i];
        let mut best_score = 0.;
        let mut last_tied = false;

        if varmer_d.last_position == 29 && varmer_d.first_position == 1 {
            let hit = &hits[i];
            log::trace!(
                "TEST HIT HIT: {:?}, score {}, bad {}",
                hit.varmer.seq_vec,
                hit.same.len(),
                hit.r_to_a.len() + hit.a_to_r.len() + hit.del.len()
            );
        }

        for j in 0..i {
            let mut score = 0;
            let mut bad = 0;
            let mut inside = false;
            for varmer in graph[&hits[i].varmer.seq_vec].in_varmers.iter() {
                if &**varmer == &hits[j].varmer.seq_vec {
                    inside = true;
                    break;
                }
            }

            if !inside {
                continue;
            }

            for pos in hits[i].same.iter() {
                if !same_vecs[j].contains(pos) {
                    score += 1;
                }
            }
            for pos in hits[i].r_to_a.iter() {
                if !rtoa_vecs[j].contains(pos) {
                    bad += 1;
                }
            }
            for pos in hits[i].a_to_r.iter() {
                if !ator_vecs[j].contains(pos) {
                    bad += 1;
                }
            }
            for pos in hits[i].del.iter() {
                if !del_vecs[j].contains(pos) {
                    bad += 1;
                }
            }

            if varmer_d.last_position == 29 && varmer_d.first_position == 1 {
                log::trace!(
                    "TEST HIT HIT j: {:?} best_index {} new scores {:?}, bad {}, score {}, INSIDE {}",
                    hits[j].varmer.seq_vec,
                    best_index,
                    best_score,
                    bad,
                    score,
                    inside
                );
            }

            if bad + rtoa_vecs[j].len() + ator_vecs[j].len() + del_vecs[j].len() <= threshold {
                let put_score = score as f64 + mismatch_pen * bad as f64 + dp_vec[j].0;
                if best_score < put_score {
                    best_index = j;
                    best_score = put_score;
                    last_tied = false;
                }
                else if best_score == put_score && !ambiguous_allowed{
                    best_index = i;
                    last_tied = true;
                }
            }
        }
        if varmer_d.last_position == 29 && varmer_d.first_position == 1 {
            log::trace!("best_index {}, i {}", best_index, i);
        }

        let j = best_index;
        if j == i || last_tied{
            continue;
        }
        traceback_vec[i] = j;
        rtoa_vecs[i] = rtoa_vecs[i].union(&rtoa_vecs[j]).cloned().collect();
        ator_vecs[i] = ator_vecs[i].union(&ator_vecs[j]).cloned().collect();
        del_vecs[i] = del_vecs[i].union(&del_vecs[j]).cloned().collect();
        same_vecs[i] = same_vecs[i].union(&same_vecs[j]).cloned().collect();
        //            ator_vecs[i] = ator_vecs[i].union(&hits[j].3).cloned().collect();
        //            del_vecs[i] = del_vecs[i].union(&hits[j].4).cloned().collect();
        //            same_vecs[i] = same_vecs[i].union(&hits[j].1).cloned().collect();
        dp_vec[i] = (best_score, dp_vec[i].1 + dp_vec[j].1);
    }

    //get traceback with max score
    //sort the vector and remember the index
    let mut score_vec = dp_vec.iter().enumerate().map(|(i, x)| (x.0, x.1, i)).collect::<Vec<(f64,u64,usize)>>();
    score_vec.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let max_index = score_vec[0].2;
    let max_score = (score_vec[0].0, score_vec[0].1);

    //When mapping reads, don't allow ambiguous paths. 
    if !ambiguous_allowed{
        if score_vec.len() > 1 && score_vec[1].0 == max_score.0{
            return empty_result;
        }
    }
    

    return DpResult {
        score: max_score,
        traceback: traceback_vec,
        max_index,
        hits,
        dels_max: del_vecs[max_index].clone(),
        rtoa_max: rtoa_vecs[max_index].clone(),
        ator_max: ator_vecs[max_index].clone(),
        same_max: same_vecs[max_index].clone(),
        total_errs: del_vecs[max_index].len()
            + rtoa_vecs[max_index].len()
            + ator_vecs[max_index].len(),
    };
}

fn fmt_varmer_d(varmer_d: &DictFrag) -> String {
    let mut seq1 = varmer_d.seq.iter().collect::<Vec<(&u32, &u8)>>();
    seq1.sort_by(|a, b| a.0.cmp(&b.0));
    let mut seq1 = seq1
        .iter()
        .map(|(_, geno)| format!("{}", *geno))
        .collect::<Vec<String>>()
        .join("");
    seq1.push_str(&format!(
        "_{}-{}--{}-LEN{}",
        varmer_d.first_position,
        varmer_d.last_position,
        varmer_d.cov,
        varmer_d.seq.len()
    ));
    seq1
}

fn print_varmer_d(varmer_d: &DictFrag) {
    let seq1 = fmt_varmer_d(varmer_d);

    log::trace!("{}", seq1);
}

pub fn query_unitigs(unitigs: &FxHashMap<VarMer, DBGInfo>, threshold: usize) -> Vec<VarMer> {
    let mut dict_unitigs = vec![];
    let mut bad_unitigs = vec![];
    for unitig in unitigs.iter() {
        let mut dict_seq = FxHashMap::default();
        let cov = unitig.1.coverage;
        let mut first_position = u32::MAX;
        let mut last_position = 0;
        for (pos, geno) in unitig.0.iter() {
            dict_seq.insert(*pos, *geno);
            if *pos < first_position {
                first_position = *pos;
            }
            if *pos > last_position {
                last_position = *pos;
            }
        }
        dict_unitigs.push(DictFrag {
            seq: dict_seq,
            seq_vec: unitig.0.clone(),
            first_position,
            last_position,
            cov,
        });
    }

    for unitig1 in dict_unitigs.iter() {
        let hits = get_hits(unitig1, &dict_unitigs, threshold, false);

        log::trace!("QUERY");
        print_varmer_d(unitig1);

        if hits.len() == 0 {
            continue;
        }

        let dp_res = dp_hits(&hits, unitig1, unitigs, threshold, 0.2, true);

        if dp_res.total_errs > threshold {
            log::trace!("FAILED ERRS: {}", dp_res.total_errs);
            continue;
        }

        if dp_res.total_errs + dp_res.same_max.len() < unitig1.seq.len() as usize {
            log::trace!("FAILED LEN: {}", dp_res.total_errs + dp_res.same_max.len());
            log::trace!(
                "VECS {:?} {:?} {:?} {:?}",
                dp_res.rtoa_max,
                dp_res.ator_max,
                dp_res.dels_max,
                dp_res.same_max
            );
            log::trace!("BEST HIT: {:?}", hits[dp_res.max_index].varmer.seq_vec);
            log::trace!("ITSELF {:?}", unitig1.seq_vec);
            continue;
        }

        log::trace!("PASSED");
        print_varmer_d(unitig1);
        log::trace!(
            "MAX SCORE: {:?}, RTOA {:?}, ATOR {:?}, DEL {:?}, max_index {}",
            dp_res.score,
            dp_res.rtoa_max,
            dp_res.ator_max,
            dp_res.dels_max,
            dp_res.max_index
        );

        log::trace!("HIT");
        print_varmer_d(hits[dp_res.max_index].varmer);
        let max_index = dp_res.max_index;
        let traceback_vec = &dp_res.traceback;

        let mut avg_cov = hits[max_index].varmer.cov;
        let mut count = 1;
        let mut final_cov = 0;

        if hits[max_index].varmer.first_position <= unitig1.first_position
            && hits[max_index].varmer.last_position >= unitig1.last_position
        {
            final_cov = hits[max_index].varmer.cov;
        }

        let mut trace_index = max_index;

        while trace_index != traceback_vec[trace_index] {
            trace_index = traceback_vec[trace_index];
            log::trace!("HIT");
            print_varmer_d(hits[trace_index].varmer);
            avg_cov += hits[trace_index].varmer.cov;
            count += 1;

            if hits[trace_index].varmer.first_position <= unitig1.first_position
                && hits[trace_index].varmer.last_position >= unitig1.last_position
            {
                final_cov = hits[max_index].varmer.cov;
            }
        }

        if final_cov == 0 {
            final_cov = avg_cov / count;
        }
        log::trace!("FINAL COV: {}, contig_COV {}", final_cov, unitig1.cov);

        if dp_res.dels_max.len() > 0 {
            if binomial_test(
                final_cov as u64,
                unitig1.cov as u64,
                0.15_f64.powi(threshold as i32),
            ) > 0.005
            {
                bad_unitigs.push(unitig1.seq_vec.clone());
                log::debug!("BAD DEL: {}", final_cov);
                print_varmer_d(unitig1);
                print_varmer_d(hits[dp_res.max_index].varmer);
            }
        }
        if dp_res.rtoa_max.len() > 0 {
            if binomial_test(
                final_cov as u64,
                unitig1.cov as u64,
                0.10_f64.powi(threshold as i32),
            ) > 0.005
            {
                bad_unitigs.push(unitig1.seq_vec.clone());
                log::debug!("BAD RTOA: {}", final_cov);
                print_varmer_d(unitig1);
                print_varmer_d(hits[dp_res.max_index].varmer);
            }
        }
        if dp_res.ator_max.len() > 0 {
            if binomial_test(
                final_cov as u64,
                unitig1.cov as u64,
                0.10_f64.powi(threshold as i32),
            ) > 0.005
            {
                bad_unitigs.push(unitig1.seq_vec.clone());
                log::debug!("BAD ATOR: {}", final_cov);
                print_varmer_d(unitig1);
                print_varmer_d(hits[dp_res.max_index].varmer);
            }
        }
    }
    return bad_unitigs;
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct Hit<'a> {
    varmer: &'a DictFrag,
    same: FxHashSet<u32>,
    r_to_a: FxHashSet<u32>,
    a_to_r: FxHashSet<u32>,
    del: FxHashSet<u32>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct DictFrag {
    seq: FxHashMap<u32, u8>,
    seq_vec: VarMer,
    first_position: u32,
    last_position: u32,
    cov: u64,
}

impl Hash for DictFrag {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seq_vec.hash(state);
    }
}

fn fragdbg_to_dictfrag(frag: &FragDBG) -> DictFrag {
    let mut seq = FxHashMap::default();
    for (pos, geno) in frag.seq.iter() {
        seq.insert(*pos, *geno);
    }
    DictFrag {
        seq,
        seq_vec: frag.seq.iter().cloned().collect(),
        first_position: frag.first_position,
        last_position: frag.last_position,
        cov: frag.seq.len() as u64,
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct AlignResult<'a> {
    same: i32,
    rtoa: i32,
    ator: i32,
    del: i32,
    path: Vec<&'a VarMer>,
    cov: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq)]
struct DpResult<'a> {
    score: (f64, u64),
    hits: &'a Vec<Hit<'a>>,
    traceback: Vec<usize>,
    max_index: usize,
    dels_max: FxHashSet<u32>,
    rtoa_max: FxHashSet<u32>,
    ator_max: FxHashSet<u32>,
    same_max: FxHashSet<u32>,
    total_errs: usize,
}

fn binomial_test(n: u64, k: u64, p: f64) -> f64 {
    // n: number of trials
    // k: number of successes
    // p: probability of success

    // Create a binomial distribution
    let binomial = Binomial::new(p, n).unwrap();

    // Calculate the probability of observing k or more successes
    let p_value = 1.0 - binomial.cdf(k);

    p_value
}

fn get_path_haps(dbg_frags: &Vec<FragDBG>, final_unitigs: &FxHashMap<VarMer, DBGInfo>, paths: Vec<Vec<(&VarMer, usize)>>) {
    let conservative = false;
    let mut path_frags: Vec<DictFrag> = vec![];
    for path in paths {
        let mut seq = vec![];
        let mut covs = vec![];
        for (varmer, cov) in path.iter() {
            for (pos, geno) in varmer.iter() {
                if !seq.contains(&(*pos, *geno)) {
                    seq.push((*pos, *geno));
                }
            }
            covs.push(*cov);
        }
        let final_cov = covs.iter().fold(0, |acc, cov| acc + cov) / covs.len();
        let mut s = seq
            .iter()
            .map(|(_, geno)| format!("{}", geno))
            .collect::<Vec<String>>()
            .join("");
        s.push_str(&format!(
            "_{}-{}:{}",
            seq.first().unwrap().0,
            seq.last().unwrap().0,
            final_cov
        ));
        log::info!("PATH: {}", s);
        let path_frag = DictFrag {
            seq: seq.iter().cloned().collect::<FxHashMap<u32, u8>>(),
            first_position: seq.first().unwrap().0,
            last_position: seq.last().unwrap().0,
            seq_vec: seq,
            cov: 0,
        };
        path_frags.push(path_frag);
    }
    let mut assignments = vec![vec![]; path_frags.len()];

    for frag in dbg_frags {
        let mut best_score = 0;
        let mut best_index = 0;
        for (i, path_frag) in path_frags.iter().enumerate() {
            let mut score = 0;
            for (pos, geno) in frag.seq.iter() {
                if path_frag.seq.contains_key(pos) {
                    if *geno == *path_frag.seq.get(pos).unwrap() {
                        score += 1;
                    }
                } else {
                    score -= 1;
                }
            }
            if score > best_score {
                best_score = score;
                best_index = i;
            } else if score == best_score {
                if path_frag.seq.len() > path_frags[best_index].seq.len() {
                    best_index = i;
                }
            }
        }
        assignments[best_index].push(frag);
    }

    //print relative percentage
    let total_cov = assignments.iter().fold(0, |acc, assignment| {
        acc + assignment.iter().fold(0, |acc, frag| acc + frag.seq.len())
    });
    for i in 0..assignments.len() {
        let mut relative_cov = 0;
        for frag in assignments[i].iter() {
            relative_cov += frag.seq.len();
        }
        log::info!(
            "PATH: {}, COV: {}, RELATIVE COV: {}",
            i,
            relative_cov,
            relative_cov as f64 / total_cov as f64
        );
    }

    //print ids to a file where each row is a path and each column is a frag id, tab sep
    let mut id_file = BufWriter::new(std::fs::File::create("ids.txt").unwrap());
    for assignment in assignments.iter() {
        for frag in assignment.iter() {
            id_file.write_all(frag.id.as_bytes()).unwrap();
            id_file.write_all(b"\t").unwrap();
        }
        id_file.write_all(b"\n").unwrap();
    }

    //Write these seq strings to a fasta file, diff identifier
    for (i, haps) in assignments.iter().enumerate() {
        let mut fasta_file =
            BufWriter::new(std::fs::File::create(format!("path_{}.fasta", i)).unwrap());
        let seqs = haps
            .iter()
            .map(|frag| frag.seq_string[0].to_ascii_vec())
            .collect::<Vec<Vec<u8>>>();
        for (j, seq) in seqs.iter().enumerate() {
            if j != 0 {
                fasta_file.write_all(b"\n").unwrap();
            }
            let rec_str = format!(">path_{}_{}\n", i, j);
            fasta_file.write_all(rec_str.as_bytes()).unwrap();
            fasta_file.write_all(seq).unwrap();
        }
    }
}


fn clean_hanging_kmers (unitigs: FxHashMap<VarMer, DBGInfo>, k: usize) -> FxHashMap<VarMer, DBGInfo> {
    let mut int_rep: FxHashMap<usize, (Vec<usize>, Vec<usize>, u64)> = FxHashMap::default();
    let unitigs_to_ints : FxHashMap<&VarMer, usize> = unitigs.iter().enumerate().map(|(i, (varmer, _))| (varmer, i)).collect();

    for (varmer, info) in unitigs.iter(){
        let mut in_ints = vec![];
        let mut out_ints = vec![];
        for in_varmer in info.in_varmers.iter(){
            in_ints.push(*unitigs_to_ints.get(&**in_varmer).unwrap());
        }
        for out_varmer in info.out_varmers.iter(){
            out_ints.push(*unitigs_to_ints.get(&**out_varmer).unwrap());
        }
        int_rep.insert(unitigs_to_ints.get(varmer).unwrap().clone(), (in_ints, out_ints, info.coverage));
    }

    let mut new_unitigs_map = FxHashMap::default();

    let mut new_unitigs = FxHashMap::default();
    let k_l = k / 2;
    let k_r = k - k_l;

    for (varmer, info) in unitigs.iter() {
        let ind = unitigs_to_ints.get(varmer).unwrap();
        let cut_varmer;
        if info.in_varmers.len() == 0 && info.out_varmers.len() == 0 {
            cut_varmer = varmer.clone();
        }
        else if info.in_varmers.len() != 0 && info.out_varmers.len() != 0 {
            cut_varmer = varmer[k_l..varmer.len() - k_r].iter().cloned().collect::<VarMer>();
        }
        else if info.in_varmers.len() == 0 && info.out_varmers.len() != 0{
            cut_varmer = varmer[..varmer.len() - k_r].iter().cloned().collect::<VarMer>();
        }
        else if info.in_varmers.len() != 0 && info.out_varmers.len() == 0{
            cut_varmer = varmer[k_l..].iter().cloned().collect::<VarMer>();
        }
        else{
            cut_varmer = varmer.clone();
        }
        new_unitigs_map.insert(ind, cut_varmer.clone());
        new_unitigs.insert(cut_varmer, DBGInfo::default());
    }

    for (ind, (in_ints, out_ints, cov)) in int_rep.iter(){
        let mut new_in_ints = vec![];
        let mut new_out_ints = vec![];
        for in_int in in_ints.iter(){
            if new_unitigs_map.contains_key(in_int){
                new_in_ints.push(*in_int);
            }
        }
        for out_int in out_ints.iter(){
            if new_unitigs_map.contains_key(out_int){
                new_out_ints.push(*out_int);
            }
        }
        let varmer = new_unitigs_map.get(ind).unwrap().clone();
        new_unitigs.get_mut(&varmer).unwrap().in_varmers = new_in_ints.iter().map(|x| Rc::new(new_unitigs_map.get(x).unwrap().clone())).collect();
        new_unitigs.get_mut(&varmer).unwrap().out_varmers = new_out_ints.iter().map(|x| Rc::new(new_unitigs_map.get(x).unwrap().clone())).collect();
        new_unitigs.get_mut(&varmer).unwrap().coverage = *cov;
    }

    return new_unitigs;
}
