use clap::Parser;
//use evalexpr::*;
//use schemars::{schema_for, JsonSchema}; 
//use serde::Deserialize;
// Note can also precompile expressions with evalexpr if desired see https://docs.rs/evalexpr/latest/evalexpr/

//see https://doc.rust-lang.org/book/ch14-02-publishing-to-crates-io.html#making-useful-documentation-comments

/// The CLI for the CTMC life cycle modelling project.
// '''
// println!("Let's model that life cycle!");
// '''

// Author: Michael P. Mariani PhD
// Year: 2025
// License: GnuPLv3
// Collaborator(s): Professor Olvia (Prosper) Feldman, PhD.,
// Associate Professor of Mathematics, University of Tennessee.

// # Some Config
// 
// This config does some stuff.
//
// # Example
//
// ```jsonc
// {
//   // use this config setting in a certain way
//   "foo": ["bar", "baz"],
// }
// ```
//#[derive(Debug, Deserialize, Default, JsonSchema)]
//pub struct Config {
//    foo: Vec<String>
//}

mod run_main;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args_help {

    /// N
    #[arg(short, long, default_value_t = 10)]
    large_n: u8,
	
	#[arg(short, long, default_value_t = 24.00/(20.00/60.00) )]
    a: f64,
	
	#[arg(short, long, default_value_t = 24.00/(20.00/60.00) )]
    b: f64,
	
	#[arg(short, long, default_value_t = 0.08)]
    r: f64,
	
	#[arg(short, long, default_value_t = 1)]
    z_mu: u8,
	
	#[arg(short, long, default_value_t = (24/19) as f64 )]
    isigma_z: f64,
	
	#[arg(short, long, default_value_t = 0.6 )]
    jsigma_e: f64,
	
	#[arg(short, long, default_value_t = 1.4 )]
    ee_mu: f64,
	
	#[arg(short, long, default_value_t = 0 )]
    o_mu: u8,
	
	#[arg(short, long, default_value_t = 3e3)]
    n_small: f64,
	
	#[arg(short, long, default_value_t = 0.2)]
    p: f64,
	
	#[arg(short, long, default_value_t = 0.39)]
    u_alpha: f64,
	
	#[arg(short, long, default_value_t = 1)]
    x_beta: u8,
	
	#[arg(short, long, default_value_t = 1)]
    y_eta: u8,
	
	#[arg(short, long, default_value_t = 8)]
    w_rho: u8,
	
	#[arg(short, long, default_value_t = (1/7) as f64 )]
    k: f64,
	
	#[arg(short, long, default_value_t = 10)]
    time_zero: u8,
	
	#[arg(short, long, default_value_t = 0.1)]
    max_bias: f64,
	
	#[arg(short, long, default_value_t = 1)]
    sim_num: u8,
	
	//par.N        = 2 				
	//par.a        = 24/(20/60) 		
	//par.b        = 24/(25/60) 		
	//par.r        = 0.08 			
	//par.mu_z     = 1 				
	//par.sigma_z  = 24/19 			
	//par.sigma_e  = 0.6 				
	//par.mu_e     = 1.4 				
	//par.mu_o     = 0 				
	//par.n        = 3e3 				
	//par.p        = 0.2 				
	//par.alpha    = 0.39 			
	//par.beta     = 1 //%0.96    	
	//par.eta      = 1 				
	//par.rho      = 8 				
	//par.k        = 1/7
	//par.t0       = 10
	//par.max_bias = 0.1 //%0.5 //%0.1
	//par.NumSim   = 1 //%10000
	
	//let matches = App::new("clap")
    //.arg(Arg::with_name("Fmin")
    //   .required(false)
    //    .takes_value(true)
    //    .short("Fmin")
    //    .multiple(false)
    //    .possible_values(&["min"])
    //)
    //.get_matches();
	
}

fn main() {
    
	let args = Args::parse();

    for _ in 0..args.sim_num {
		
        println!("Number of simulations {}!", args.sim_num);
		
	}
	
	let args: Vec<String> = env::args().collect();
    dbg!(args);
	
	
}

//struct Cli {
//    // Define a command field to hold the command to execute
//    #[clap(name = "command")]
//    command: String,
//
//    // Define an args field to hold arguments for the command
//    #[clap(name = "args")]
//    args: Vec<String>,
//
//    // Define an append field to specify whether to append to an output file
//    #[clap(
//        short = 'a',
//        long,
//        help = "Append to output file instead of overwriting"
//    )]
//    append: bool,
//
//    // Define an output field to specify an output file
//    #[clap(short = 'o', long, help = "Write output to a file instead of stdout")]
//    output: Option<String>,
//}

//cargo run -- sleep 2
//cargo run -- -o output.txt sleep 10
