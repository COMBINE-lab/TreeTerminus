use petgraph::unionfind::UnionFind;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::fs::*;
use std::io::Write;
use std::io::{self, BufReader};
use std::path::PathBuf;
extern crate serde_json;
extern crate serde_pickle;
extern crate serde_stacker;

use std::ffi::CString;
use std::os::raw::{c_char, c_int};

use crate::binary_tree::{get_binary_rooted_newick_string, sort_group_id, TreeNode};
use crate::salmon_types::{ConsensusFileList, FileList};

extern "C" {
    pub fn run_cons(argsc:c_int, argsv:*const *const c_char) -> c_int;
}

#[allow(unused_variables)]
fn create_union_find(g: &[String], ntxps: usize) -> UnionFind<usize> {
    let mut unionfind_struct = UnionFind::new(ntxps);
    let mut visited: Vec<i32> = vec![-1; ntxps];
    let mut count = 0;
    for (_i, group) in g.iter().enumerate() {
        let g_set: Vec<usize> = group
            .clone()
            .split('_')
            .map(|x| x.parse::<usize>().unwrap())
            .collect();
        let source = g_set[0];
        if visited[source as usize] == -1 {
            count += 1;
            visited[source as usize] = 0
        }
        for t in g_set.iter().skip(1) {
            unionfind_struct.union(source, *t);
            if visited[*t as usize] == -1 {
                count += 1;
                visited[*t as usize] = 0
            }
        }
    }
    // if count != ntxps {
    //     panic!("The number of expected transcripts {} do not match the counted transcripts from the groups{}", ntxps, count);
    // }
    unionfind_struct
}

fn get_merged_bparts(
    groups: &HashMap<usize, Vec<usize>>,
    all_groups_bpart: &HashMap<String, HashMap<String, u32>>,
    uf: &UnionFind<usize>,
) -> HashMap<String, HashMap<String, u32>> {
    let all_groups: Vec<String> = all_groups_bpart.keys().cloned().collect();
    let mut merged_bparts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    for (_j, old_g) in all_groups.iter().enumerate() {
        let f_txp = old_g
            .clone()
            .split('_')
            .map(|x| x.parse::<usize>().unwrap())
            .collect::<Vec<usize>>()[0];
        let g_ind = uf.find(f_txp); //the vertex index
        let strings: Vec<String> = groups[&g_ind].iter().map(|n| n.to_string()).collect();
        let m_group = strings.join("_").to_string();
        let m_bpart_key = merged_bparts
            .entry(sort_group_id(&m_group.clone()))
            .or_insert_with(HashMap::new);

        for (b_part, count) in all_groups_bpart.get(&old_g.clone()).unwrap().iter() {
            let c_count = m_bpart_key.entry(b_part.clone()).or_insert(0);
            *c_count += count;
        }
    }
    merged_bparts
}

#[allow(dead_code)]
fn find_groups_in_merged(
    groups: &HashMap<usize, Vec<usize>>,
    all_groups: &[String],
    uf: &UnionFind<usize>,
) -> HashMap<String, Vec<String>> {
    // Given a grouped hashmap that contains all child txps under a main txp (together will form a (merged) group),
    // vector of all keys of a group and a union find, returns a hashmap containing merged group and old groups
    let mut merged_groups: HashMap<String, Vec<String>> = HashMap::new();
    for (j, old_g) in all_groups.iter().enumerate() {
        let f_txp = old_g
            .clone()
            .split('_')
            .map(|x| x.parse::<usize>().unwrap())
            .collect::<Vec<usize>>()[0];
        let g_ind = uf.find(f_txp); //the vertex index
        let strings: Vec<String> = groups[&g_ind].iter().map(|n| n.to_string()).collect();
        let m_group = strings.join("_").to_string();
        merged_groups
            .entry(m_group)
            .or_insert_with(Vec::new)
            .push(all_groups[j].clone());
    }
    merged_groups
}

#[allow(dead_code)]
pub fn merge_groups(
    all_groups_bpart: &HashMap<String, HashMap<String, u32>>,
    ntxps: usize,
) -> HashMap<String, HashMap<String, u32>> {
    let all_groups: Vec<String> = all_groups_bpart.keys().cloned().collect();
    let g_union = create_union_find(&all_groups, ntxps as usize);
    let mut groups = HashMap::new();
    for i in 0..ntxps {
        let root = g_union.find(i);
        if root != i {
            groups.entry(root).or_insert_with(|| vec![root]).push(i);
        }
    }
    get_merged_bparts(&groups, all_groups_bpart, &g_union)
}

fn comp_diff(m_group: &str, oth_groups: &[String]) -> Vec<usize> {
    let par_set: HashSet<usize> = m_group
        .split('_')
        .map(|x| x.parse::<usize>().unwrap())
        .collect();
    let mut child_set: HashSet<usize> = HashSet::new();
    for g in oth_groups.iter() {
        let c: HashSet<usize> = g.split('_').map(|x| x.parse::<usize>().unwrap()).collect();
        child_set.extend(&c);
    }
    let diff: Vec<usize> = par_set.difference(&child_set).cloned().collect();
    diff
}

fn get_group_trees(
    merged_group: &String,
    groups: &[String],
    samp_group_trees: &[HashMap<String, TreeNode>],
) -> (String, Vec<String>) {
    // Returns a tuple that contains merged group and total sample count along with child group and the number of sample that child group appears in
    let mut samp_nwk: Vec<String> = Vec::new();
    let mut g_inf = String::from("");
    for (_i, samp_hash) in samp_group_trees.iter().enumerate() {
        let mut g_vec: Vec<String> = Vec::new();
        let mut s_trees: Vec<String> = Vec::new();
        for (_j, g) in groups.iter().enumerate() {
            if samp_hash.contains_key(g) {
                g_vec.push(g.clone());
                //println!("{}\t{:?}",g, samp_group_trees[_i].get(g).unwrap().traverse_tree());
                s_trees.push(get_binary_rooted_newick_string(
                    samp_group_trees[_i].get(g).unwrap(),
                ));
            }
        }
        let l = s_trees.len();
        let mut cur_nwk_trees = match l {
            0 => "(".to_string(),
            1 => s_trees[0].clone(),
            _ => format!("({}", s_trees.join(",")),
        };

        // let mut cur_nwk_trees = "".to_string();
        // if length(s_trees) <= 1 {
        //     cur_nwk_trees.insert_str(0, "(");
        //     if length(s_trees) == 1 {
        //         cur_nwk_trees.push_str(&format!("{},",s_trees[0]));
        //     }
        // }
        // let mut cur_nwk_trees = format!("({}",s_trees.join(","));
        let diff = comp_diff(merged_group, &g_vec);
        if !diff.is_empty() {
            if l == 1 {
                cur_nwk_trees.insert(0, '(');
            }
            if l > 0 {
                for v in diff.iter() {
                    cur_nwk_trees.push_str(&format!(",{}", v));
                }
            } else {
                cur_nwk_trees.push_str(&format!("{}", diff[0].clone()));
                for v in diff.iter().skip(1) {
                    cur_nwk_trees.push_str(&format!(",{}", v));
                }
            }
        }

        if l != 1 || !diff.is_empty() {
            cur_nwk_trees.push_str(");");
        } else {
            cur_nwk_trees.push(';');
        }
        samp_nwk.push(cur_nwk_trees.clone());
        let gs = g_vec.join(",");
        g_inf.push_str(&format!("\t{}\t{}", _i, gs));
    }
    //    println!("{}\t{}", merged_group.clone(), count);
    g_inf.insert_str(0, &merged_group.to_string());
    (g_inf, samp_nwk)
}

fn write_file(f: &mut File, st: String) -> Result<bool, io::Error> {
    writeln!(f, "{}", st)?;
    Ok(true)
}

fn get_cons(out: &String, samp_trees: &[String]) -> String {
    let dir = PathBuf::from(out);
    let inp_nwk = dir.as_path().join("inp_tree.nwk");
    let out_nwk = dir.as_path().join("out_tree.nwk");

    let mut f_inp = File::create(inp_nwk.clone()).expect("could not create input newick file");
    for g in samp_trees.iter() {
        let _t = write_file(&mut f_inp, g.clone());
    }

    let mut args: Vec<CString> = Vec::new();
    args.push(CString::new(inp_nwk.as_path().display().to_string()).unwrap());
    args.push(CString::new(out_nwk.as_path().display().to_string()).unwrap());
    args.push(CString::new("mre").unwrap());
    args.push(CString::new("0").unwrap());
    args.push(CString::new("0").unwrap());
    args.push(CString::new("0").unwrap());
    args.push(CString::new("0").unwrap());

    let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
    let args_len: c_int = arg_ptrs.len() as c_int;
    unsafe{run_cons(args_len, arg_ptrs.as_ptr())};

    // // println!("{:?}", samp_trees);
    // let (_code, _output, _error) = run_script::run_script!(
    //     r#"
    //     ./phylip_consensus/consense < phylip_consensus/input
    //      exit 0
    //      "#
    // )
    // .unwrap();
    let cons_nwk = read_to_string(out_nwk.clone()).expect("Something went wrong reading the file");
    let (_code, _output, _error) = run_script::run_script!(
        &format!("rm {}", out_nwk.as_path().display().to_string())
    )
    .unwrap();
    // println!("{}", cons_nwk);
    // let (_code, _output, _error) = run_script::run_script!(
    //     r#"
    //     rm out*
        
    //         exit 0
    //         "#
    // )
    // .unwrap();
    cons_nwk
}

pub fn use_phylip(dir_paths: &[&str], out: &String, all_groups: &[String], ntxps: usize) {
    let g_union = create_union_find(all_groups, ntxps as usize);
    let mut groups = HashMap::new();
    for i in 0..ntxps {
        let root = g_union.find(i);
        if root != i {
            groups.entry(root).or_insert_with(|| vec![root]).push(i);
        }
    }
    let mg = find_groups_in_merged(&groups, all_groups, &g_union); //merged groups
    println!("Length of groups after merging {}", mg.len());

    println!("Reading group trees");
    let mut samp_group_trees: Vec<HashMap<String, TreeNode>> = Vec::new(); //Vector containing group trees from each sample
    let mut msamp_nwk_file: Vec<File> = Vec::new(); //Vector containing newick trees corresponding to each group
                                                    // Storing group trees in each sample in an array along with ....
    for (_i, dname) in dir_paths.iter().enumerate() {
        let compo: Vec<&str> = dname.rsplit('/').collect();
        let experiment_name = compo[0];
        let mut prefix_path = out.clone();
        prefix_path.push('/');
        prefix_path.push_str(experiment_name);

        let file_list_out = FileList::new(prefix_path);
        msamp_nwk_file.push(
            File::create(file_list_out.mgroup_nwk_file).expect("Could not open mgroup nwk file"),
        );
        let file = File::open(file_list_out.collapse_order_file);
        let reader = BufReader::new(file.unwrap());

        let mut deserializer = serde_json::Deserializer::from_reader(reader);
        deserializer.disable_recursion_limit();
        samp_group_trees.push(HashMap::deserialize(&mut deserializer).unwrap());
        //Pushing hashmap containing trees corresponding to each group
    }
    println!("Finished reading group trees");
    let file_list_out = ConsensusFileList::new(out.clone());
    let mut mg_file =
        File::create(file_list_out.merged_groups_file).expect("could not create merged group file");
    let mut clust_nwk_file =
        File::create(file_list_out.cons_nwk_file).expect("could not create cluster newick file");

    let inp_nwk_s = format!("{}/inp_tree.nwk", out.clone());
    // let (_code, _output, _error) =
    //     run_script::run_script!(&format!("echo {}  > phylip_consensus/input", inp_nwk_s)).unwrap();
    // let (_code, _output, _error) = run_script::run_script!(
    //     r#"
    //     echo "R Yes" >> phylip_consensus/input
    //     echo "Y"  >> phylip_consensus/input
    //      exit 0
    //      "#
    // )
    // .unwrap();
    for (merged_group, old_group) in mg {
        let group_inf = get_group_trees(&merged_group, &old_group, &samp_group_trees); //
        let _t = write_file(&mut mg_file, group_inf.0);
        println!("Computing cluster for group {}", merged_group.clone());
        for (_i, g) in group_inf.1.iter().enumerate() {
            let _t = write_file(&mut msamp_nwk_file[_i], g.clone());
        }
        //println!("{:?}", group_inf.1);
        //println!("{}", get_cons(out, &group_inf.1));
        let _t = write_file(&mut clust_nwk_file, get_cons(out, &group_inf.1));
    }
    let (_code, _output, _error) = run_script::run_script!(
        &format!("rm {}", inp_nwk_s)
    )
    .unwrap();
}
