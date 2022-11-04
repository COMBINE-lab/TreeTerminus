pub mod binary_tree;
mod collapse;
pub mod salmon_types;
mod util;

extern crate serde_json;
extern crate serde_pickle;
extern crate serde_stacker;

use std::collections::HashMap;
use std::fs::*;
use std::io::Write;
use std::io::{self, BufReader};

use clap::{App, AppSettings, Arg, ArgMatches, SubCommand};
use ndarray::prelude::*;
use num_format::{Locale, ToFormattedString};

use petgraph::algo::connected_components;
use petgraph::unionfind::UnionFind;
use std::path::PathBuf;

use serde::Deserialize;
use serde_json::json;

use crate::salmon_types::ConsensusFileList;

// Name of the program, to be used in diagnostic messages.
static PROGRAM_NAME: &str = "treeterminus";

/// Exit the program, printing an error message on stderr, and returning
/// a specific error code. The program name is prefixed onto the front of
/// the error message. If logging is enabled we also write the error
/// message to the log file.
#[allow(unused)]
fn exit_with_error(status: i32, message: &str) {
    //error!(target: "log_messages", "{}", message);
    writeln!(
        &mut std::io::stderr(),
        "{} ERROR: {}!",
        PROGRAM_NAME,
        message
    )
    .unwrap();
    std::process::exit(status);
}

fn do_group(sub_m: &ArgMatches) -> Result<bool, io::Error> {
    //let mut groups: Vec<Vec<usize>> = Vec::new();
    //let dir_paths : Vec<_> = sub_m.values_of("dirs").unwrap().collect();
    let dname: String = sub_m.value_of("dir").unwrap().to_string();
    let prefix: String = sub_m.value_of("out").unwrap().to_string();
    create_dir_all(prefix.clone())?;

    // let mut unionfind_vec = Vec::<UnionFind<_>>::with_capacity(dir_paths.len()) ;
    // let mut gibbs_array_vec = Vec::<Array2<f64>>::with_capacity(dir_paths.len());
    // let mut meta_info_array = Vec::<salmon_types::MetaInfo>::with_capacity(dir_paths.len());
    // let mut num_global_targrts = 0u32 ;

    let seed = sub_m
        .value_of("seed")
        .unwrap()
        .parse::<u64>()
        .expect("generate random values from the seed");
    let min_spread = sub_m
        .value_of("min-spread")
        .unwrap()
        .parse::<f64>()
        .expect("could not convert min-spread to float value");

    let tolerance = sub_m
        .value_of("tolerance")
        .unwrap()
        .parse::<f64>()
        .expect("could not convert tolerance to float value");

    let thr_bool = sub_m
        .value_of("thresh")
        .unwrap()
        .parse::<bool>()
        .expect("could not parse thr_bool");

    let mean_inf = sub_m
        .value_of("mean")
        .unwrap()
        .parse::<bool>()
        .expect("could not parse mean");

    let inf_perc = sub_m
        .value_of("inf_perc")
        .unwrap()
        .parse::<f64>()
        .expect("could not parse inf percentile");

    let mut dir_paths: Vec<String> = Vec::new();
    if mean_inf {
        let sd = read_dir(dname.clone());
        for f in sd.unwrap() {
            let x = f
                .unwrap()
                .path()
                .to_str()
                .unwrap()
                .rsplit('/')
                .collect::<Vec<_>>()[0]
                .to_string();
            if x == "quant.sf" {
                panic!("A directory above this level is required or flag mean_inf should be set to false");
            }
        }
        let sal_dir_paths = read_dir(dname.clone())?
            .map(|res| res.map(|e| e.path()))
            .filter(|res| res.as_ref().unwrap().is_dir())
            .collect::<Result<Vec<_>, io::Error>>()?;

        for entry in sal_dir_paths.iter() {
            dir_paths.push(entry.as_path().to_str().unwrap().to_string());
        }
    }

    println!("------input configuration------");
    println!("seed : {}", seed);
    println!("min-spread : {}", min_spread);
    println!("tolerance : {}", tolerance);
    println!("dir : {}", dname.clone());
    let compo: Vec<&str> = dname.rsplit('/').collect();
    //println!("{:?}",compo);
    let experiment_name = compo[0];
    let mut prefix_path = prefix;
    let file_list = salmon_types::FileList::new(dname.to_string());

    let mut file_list_vec = Vec::new();
    if !mean_inf {
        prefix_path.push('/');
        prefix_path.push_str(experiment_name);
    }

    if mean_inf {
        for dir in dir_paths {
            file_list_vec.push(salmon_types::FileList::new(dir.to_string()));
        }
    }

    // create output directory

    println!("output folder: {}", prefix_path.clone());
    println!("------------------------------");
    // create
    create_dir_all(prefix_path.clone())?;
    let file_list_out = salmon_types::FileList::new(prefix_path.clone());

    // Load the gibbs samples
    let mut x;
    let mut gibbs_array = Array2::<f64>::zeros((1, 1));
    let mut gibbs_array_vec = Vec::new();
    let mut x_vec = Vec::new();
    let mut gibbs_mat_mean = Array1::<f64>::zeros(1);
    let mut eq_class_counts: Vec<u32> = Vec::new();
    let mut eq_class;
    let mut l = 0;

    // Think about enum representation
    if mean_inf {
        let mut eq_class_vec = Vec::new();
        for (_i, file_list) in file_list_vec.iter().enumerate() {
            x_vec.push(util::parse_json(&file_list.mi_file).unwrap());
            gibbs_array_vec.push(Array2::<f64>::zeros((
                x_vec[_i].num_valid_targets as usize,
                x_vec[_i].num_bootstraps as usize,
            )));
            #[allow(unused_assignments)]
            if _i == 0 {
                gibbs_mat_mean = Array1::<f64>::zeros(gibbs_array_vec[0].shape()[0] as usize);
                x = x_vec[0].clone();
            }
            util::read_gibbs_array(
                &file_list.bootstrap_file,
                &x_vec[_i],
                &mut gibbs_array_vec[_i],
            );
            gibbs_mat_mean += &gibbs_array_vec[_i].mean_axis(Axis(1)).unwrap();
            println!("parsing eqfile {:?}", file_list.eq_file);
            eq_class_vec.push(util::parse_eq(&file_list.eq_file).unwrap());

            println!("length of eqclass {:?}", eq_class_vec[_i].neq);
            eq_class_counts.extend(vec![0_u32; eq_class_vec[_i].neq]);

            for (j, eq) in eq_class_vec[_i].classes.iter().enumerate() {
                eq_class_counts[l + j] = eq.2;
                if j == eq_class_vec[_i].neq - 1 {
                    l += j;
                }
            }
            l += 1;
        }
        x = x_vec[0].clone();
        gibbs_mat_mean /= gibbs_array_vec.len() as f64;
        eq_class = eq_class_vec[0].clone();
        #[allow(clippy::needless_range_loop)]
        for j in 1..eq_class_vec.len() {
            eq_class.neq += eq_class_vec[j].neq;
            let prev_off = eq_class.classes.offsets[eq_class.classes.offsets.len() - 1];
            let mut offset = vec![0_usize; eq_class_vec[j].classes.offsets.len() - 1];
            #[allow(clippy::needless_range_loop)]
            for _i in 0..offset.len() {
                offset[_i] = eq_class_vec[j].classes.offsets[_i + 1] + prev_off;
            }
            eq_class.classes.offsets.extend(offset.iter().copied());
            eq_class
                .classes
                .labels
                .extend(eq_class_vec[j].classes.labels.iter().copied());
            eq_class
                .classes
                .weights
                .extend(eq_class_vec[j].classes.weights.iter().copied());
            eq_class
                .classes
                .counts
                .extend(eq_class_vec[j].classes.counts.iter().copied());
        }
        // println!("{}, {}, {}, {}", eq_class.neq, eq_class.classes.offsets.len(),
        //     eq_class.classes.labels.len(), eq_class.classes.weights.len());
        // println!("{}", eq_class_counts.len());
        // println!("{}", eq_class.classes.offsets[344782]);
    } else {
        x = util::parse_json(&file_list.mi_file).unwrap();
        gibbs_array =
            Array2::<f64>::zeros((x.num_valid_targets as usize, x.num_bootstraps as usize));
        util::read_gibbs_array(&file_list.bootstrap_file, &x, &mut gibbs_array);
        gibbs_mat_mean = gibbs_array.mean_axis(Axis(1)).unwrap();

        println!("parsing eqfile {:?}", file_list.eq_file);
        eq_class = util::parse_eq(&file_list.eq_file).unwrap();
        // println!("length of eqclass {:?}", eq_class.neq);
        eq_class_counts = vec![0_u32; eq_class.neq];
        // let mut i = 0_usize;
        for (i, eq) in eq_class.classes.iter().enumerate() {
            eq_class_counts[i] = eq.2;
        }
    }

    // if a2g exists also dumps gene level groups
    let allele2txp = PathBuf::from(sub_m.value_of("a2t").unwrap().to_string());
    let asemode: bool = allele2txp.as_path().is_file();
    if asemode {
        println!(
            "Alleles would be collapsed according to the file: {:?}",
            allele2txp.to_str().unwrap()
        );
    }

    // if t2g exists restrict equivalence classes to gene level groups
    let transcript2gene = PathBuf::from(sub_m.value_of("t2g").unwrap().to_string());
    let txpmode: bool = transcript2gene.as_path().is_file();
    if txpmode {
        println!(
            "Txps within a gene would be collapsed using : {:?}",
            transcript2gene.to_str().unwrap()
        );
    }

    // take the transcript to gene mapping
    // this will also create a map from transcript id
    // to gene id

    let mut collapse_order: Vec<binary_tree::TreeNode> = Vec::new();
    for i in 0..eq_class.ntarget {
        collapse_order.push(binary_tree::TreeNode::create_leaf(i.to_string()));
    }

    // fill targets from eq_class
    let mut tnames: HashMap<String, usize> = HashMap::new();
    for i in 0..eq_class.ntarget {
        tnames.insert(eq_class.targets[i].clone(), i);
    }
    let ntarget = eq_class.ntarget;

    let mut gene2allele_map: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut allele2gene_map: Vec<usize> = Vec::new();

    let mut txp2allele_map: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut allele2txp_map: Vec<usize> = Vec::new();

    if asemode {
        util::get_map_bw_ent(
            &allele2txp,
            &mut txp2allele_map,
            &mut allele2txp_map,
            &tnames,
        );
        let nalleles = allele2txp_map.len();
        if nalleles != ntarget {
            panic!(
                "number of alleles {} not equal to number of txps in eq class file {}",
                nalleles, ntarget
            );
        }
    }

    if txpmode {
        util::get_map_bw_ent(
            &transcript2gene,
            &mut gene2allele_map,
            &mut allele2gene_map,
            &tnames,
        );
        let ntxps = allele2gene_map.len();
        if ntxps != ntarget {
            panic!(
                "number of txps {} not equal to number of txps in eq class file {}",
                ntxps, ntarget
            );
        }
    }

    if asemode {
        let nalleles = allele2txp_map.len();
        if txpmode {
            let ntxps = allele2gene_map.len();
            if nalleles != ntxps {
                panic!(
                    "number of alleles {} not equal to number of txps in eq class file {}",
                    nalleles, ntxps
                );
            }
        }
    }

    // let inf_perc = 0.25f64;
    let p = match mean_inf {
        false => util::get_infrv_percentile(&gibbs_array, inf_perc),
        true => {
            let mut cum_infrv_perc: f64 = 0.0;
            for (_i, gb) in gibbs_array_vec.iter().enumerate() {
                let perc = util::get_infrv_percentile(gb, inf_perc);
                if _i == 0 {
                    cum_infrv_perc = perc;
                } else {
                    cum_infrv_perc = cum_infrv_perc.min(perc);
                }
            }
            cum_infrv_perc
        }
    };

    println!("the {}% of infRV was : {}", inf_perc * 100., p);

    let thr = match thr_bool {
        true => {
            if !mean_inf {
                util::get_threshold(&gibbs_array, p, seed, &file_list_out)
            } else {
                let mut thresh = 0.0;
                for gb in gibbs_array_vec.iter() {
                    thresh += util::get_threshold(gb, p, seed, &file_list_out);
                }
                thresh / (gibbs_array_vec.len() as f64)
            }
        }
        false => 1e7,
    };

    println!("threshold: {}", thr);
    println!("{}", eq_class.ntarget);

    //let dpath = Path::new(file_list_out.delta_file.clone());
    let mut dfile =
        File::create(file_list_out.delta_file.clone()).expect("could not create collapse.log");
    let mut unionfind_struct: UnionFind<usize> = UnionFind::new(eq_class.ntarget);
    let mut group_order: Vec<String> = Vec::with_capacity(eq_class.ntarget);
    for i in 0..eq_class.ntarget {
        group_order.push(i.to_string())
    }

    // pass the gene to transcript mapping to the building graph phase to
    // restrict the creation of two edge between nodes from the same gene
    let mut gcfile = File::create(file_list_out.golden_collapses_log_file.clone())
        .expect("could not create golden collapse.log");
    let mut allele_file = File::create(file_list_out.allele_collapses_log_file.clone())
        .expect("could not create golden collapse.log");
    let mut gr = util::eq_experiment_to_graph(
        &eq_class,
        &mut gibbs_array,
        &mut gibbs_array_vec,
        &eq_class_counts,
        tolerance,
        thr,
        p,
        min_spread,
        &mut dfile,
        &mut unionfind_struct,
        &allele2gene_map,
        &txp2allele_map,
        &mut group_order,
        &mut collapse_order,
        mean_inf,
        &mut gcfile,
        &mut allele_file,
    );

    util::verify_graph(&eq_class_counts, &mut gr);
    // Go over the graph and keep collapsing
    // edges until we hit a point where there
    // are no more edges to that satisfies the criteria
    // and we collapse

    // connected coponents
    let num_connected_components = connected_components(&gr);
    println!("#Connected components {:?}", num_connected_components);

    let mut num_collapses = 0_usize;

    //let cpath = Path::new(file_list_out.collapsed_log_file.clone());
    let mut cfile = File::create(file_list_out.collapsed_log_file.clone())
        .expect("could not create collapse.log");

    let gcomp: Vec<petgraph::prelude::NodeIndex> = gr
        .node_indices()
        .map(|x| petgraph::graph::NodeIndex::new(x.index()))
        .collect();
    util::work_on_component(
        &eq_class_counts,
        &mut gibbs_array,
        &mut gibbs_array_vec,
        &mut gibbs_mat_mean,
        &mut unionfind_struct,
        &mut gr,
        &gcomp,
        &mut num_collapses,
        thr,
        p,
        &mut cfile,
        &mut collapse_order,
        mean_inf,
    );

    // //write down the groups
    let mut groups = HashMap::new();
    //let mut grouped_set = HashSet::new();
    for i in 0..(x.num_valid_targets as usize) {
        let root = unionfind_struct.find(i);
        if root != i {
            groups.entry(root).or_insert_with(Vec::new).push(i);
        }
    }

    // println!("Number of collapsed transcripts from conn components with 2 {}", num_collapses_2.to_formatted_string(&Locale::en));
    println!(
        "Number of collapses {}",
        num_collapses.to_formatted_string(&Locale::en)
    );

    let params = json!({
        "seed":seed,
        "tolerance":tolerance,
        "mean_inf":mean_inf,
        "thr_bool":thr_bool,
        "inp_dir":dname.clone(),
        "out_dir":prefix_path.clone(),
        "allele_mode":asemode,
        "txp_mode":txpmode,
        "inf_perc":inf_perc,
        "p":p,
        "thr":thr,
        "ntxps":eq_class.ntarget,
        "connected_components":num_connected_components,
        "ncollapses":num_collapses,
    });

    let param_log_file =
        File::create(file_list_out.param_log_file).expect("could not create group order file");
    serde_json::to_writer(param_log_file, &params)?;
    let mut gfile = File::create(file_list_out.group_file).expect("could not create groups.txt");
    let mut co_file = File::create(file_list_out.collapse_order_file)
        .expect("could not create collapse order file");
    let nwk_path = match mean_inf {
        true => {
            let cons = ConsensusFileList::new(prefix_path.clone());
            cons.cons_nwk_file
        }
        false => file_list_out.group_nwk_file,
    };
    let mut nwk_file = File::create(nwk_path).expect("could not create group write file");
    let _write = util::group_writer(&mut gfile, &groups);
    let _write = util::collapse_order_writer(&mut co_file, &mut nwk_file, &groups, &collapse_order);

    Ok(true)
}

fn do_collapse(sub_m: &ArgMatches) -> Result<bool, io::Error> {
    let sal_dir: String = sub_m.value_of("dirs").unwrap().to_string();
    let _md = metadata(sal_dir.clone()).unwrap_or_else(|_| panic!("Invalid directory {}", sal_dir));

    let sal_dir_paths = read_dir(sal_dir)?
        .map(|res| res.map(|e| e.path()))
        .filter(|res| res.as_ref().unwrap().is_dir())
        .collect::<Result<Vec<_>, io::Error>>()?;

    let mut dir_paths: Vec<&str> = Vec::new();
    for entry in sal_dir_paths.iter() {
        dir_paths.push(entry.as_path().to_str().unwrap());
    }
    let prefix: String = sub_m.value_of("out").unwrap().to_string();

    let mut bipart_counter: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut group_keys: Vec<String> = Vec::new();
    // let mut l = 0; // num of groups in 1st
    let mut ntxps = 0; // num of transcripts
    let mut tnames: Vec<String> = Vec::new();
    // add edges
    for (i, dname) in dir_paths.iter().enumerate() {
        //let dname = dname.as_path().to_str().unwrap();
        let compo: Vec<&str> = dname.rsplit('/').collect();
        if compo[0].starts_with('.') {
            continue;
        }
        let experiment_name = compo[0];
        println!("experiment name {}", experiment_name);

        let mut dir_bipart_counter: HashMap<String, HashMap<String, u32>> = HashMap::new(); // Storing counts of each bipartition

        let mut prefix_path = prefix.clone();
        prefix_path.push('/');
        prefix_path.push_str(experiment_name);
        create_dir_all(prefix_path.clone())?;
        let file_list_out = salmon_types::FileList::new(prefix_path);
        let file = File::open(file_list_out.collapse_order_file);
        let reader = BufReader::new(file.unwrap());

        let mut deserializer = serde_json::Deserializer::from_reader(reader);
        deserializer.disable_recursion_limit();

        //println!("{:?}", deserializer);
        let collapse_order: HashMap<String, binary_tree::TreeNode> =
            HashMap::deserialize(&mut deserializer).unwrap(); // can be replaced by vector
        for key in collapse_order.keys() {
            let node = collapse_order.get(key).unwrap();
            let req_group = binary_tree::sort_group_id(&node.id);
            //let node_vec = group_bipart.entry(node.id.clone()).or_insert(Vec::<String>::new());
            let dir_group_key = dir_bipart_counter
                .entry(req_group.clone())
                .or_insert_with(HashMap::new);
            let overall_group_key = bipart_counter
                .entry(req_group.clone())
                .or_insert_with(HashMap::new);

            //binary_tree::compute_bipart_count(node, &mut bipart_counter, &mut dir_bipart_counter, &node_set, node_vec);
            group_keys.push(req_group.clone());
            binary_tree::compute_bipart_count2(node, overall_group_key, dir_group_key);
        }
        if i == 0 {
            let file_old = salmon_types::FileList::new(dname.to_string());
            let eq_class = util::parse_eq(&file_old.eq_file).unwrap();
            ntxps = eq_class.ntarget;

            // for i in 0..eq_class.targets.len()
            // {
            //     tnames.push(i.to_string());
            // }
            tnames.extend(eq_class.targets.clone());
        }
        println!(
            "Number of groups in {} are {}",
            dname,
            dir_bipart_counter.len()
        );
        let mut bipart_file = File::create(file_list_out.group_bp_splits_file)
            .expect("could not create group bp splits");
        let _f = util::MapTrait::bipart_writer(&dir_bipart_counter, &mut bipart_file, &tnames);
    }
    let all_groups: Vec<String> = bipart_counter.keys().cloned().collect();
    collapse::use_phylip(&dir_paths, &prefix, &all_groups, ntxps);

    // filter based on the threshold

    Ok(true)
}

// The entry point of the program
// that reads from the files and build
// graphs or later produces the collapsed
// files.

fn main() -> io::Result<()> {
    let matches = App::new("TreeTerminus")
	.setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1.0")
        .author("Singh et al.")
        // .about("Data-driven grouping of transcripts to reduce inferential uncertainty")
        .subcommand(
            SubCommand::with_name("group")
            .about("perform per-sample grouping of transcripts; required prior to consensus collapse.")
            .arg(
                Arg::with_name("dir")
                    .long("dir")
                    .short("d")
                    .required(true)
                    .takes_value(true)
                    .help("directory to read input from")
            )
            .arg(
                Arg::with_name("min-spread")
                    .long("min-spread")
                    .short("m")
                    .takes_value(true)
                    .default_value("0.1")
                    .help("the minimum spread a transcript must exhibit to enable an \
                          attached edge to be a collapse candidate")
            )
            .arg(
                Arg::with_name("tolerance")
                    .long("tolerance")
                    .takes_value(true)
                    .default_value("0.001")
                    .help("The allowable difference between the weights of transcripts \
                          in same equivalence classes to treat them as identical")
            )
            .arg(
                Arg::with_name("seed")
                    .long("seed")
                    .takes_value(true)
                    .default_value("10")
                    .help("seed for random generator")
            )
            .arg(
                Arg::with_name("a2t")
                    .long("a2t")
                    .takes_value(true)
                    .default_value("")
                    .help("Mapping allele to transcript")
            )
            .arg(
                Arg::with_name("t2g")
                    .long("t2g")
                    .takes_value(true)
                    .default_value("")
                    .help("Mapping transcript to gene")
            )
            .arg(
                Arg::with_name("thresh")
                    .long("thr")
                    .takes_value(true)
                    .default_value("false")
                    .help("should a threshold be used for collapsing")
            )
            .arg(
                Arg::with_name("out")
                    .long("out")
                    .short("o")
                    .required(true)
                    .takes_value(true)
                    .requires("dir")
                    .help("prefix where output would be written")
            )
            .arg(
                Arg::with_name("mean")
                    .long("mean_inf")
                    .takes_value(true)
                    .default_value("true")
                    .help("mean infrv for tree construction")
            )
            .arg(
                Arg::with_name("inf_perc")
                .long("inf_perc")
                .short("i")
                .takes_value(true)
                .default_value("0")
                .help("inferential variance percentile threshold that determines whether a transcript will be considered for grouping")
            )
        )
        .subcommand(
            SubCommand::with_name("consensus")
            .about("Produce a set of consensus trees from the individual per-sample trees obtained for an RNA-Seq experiment after running the group step.")
            .arg(
                Arg::with_name("dirs")
                    .long("dirs")
                    .short("d")
                    .required(true)
               //     .multiple(true)
                    .takes_value(true)
                    .help("direcotories to read the group files from")
            )
            .arg(
                Arg::with_name("out")
                    .long("out")
                    .short("o")
                    .required(true)
                    .takes_value(true)
                    .requires("dirs")
                    .help("prefix where output would be written")
            )
        ).get_matches();

    pretty_env_logger::init_timed();

    match matches.subcommand() {
        ("group", Some(sub_m)) => {
            do_group(sub_m).expect("Grouping failed");
        }
        ("consensus", Some(sub_m)) => {
            do_collapse(sub_m).expect("Grouping failed");
        }
        _ => unreachable!(),
    }

    Ok(())
}
