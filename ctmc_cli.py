# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 17:01:25 2025

@author: Michael P. Mariani PhD
"""

import json, sys, click

@click.command(help="Python CLI for the CTMC Life Cycle Project")
#@click.option("-i", "--input", "infile", type=click.File(), default=sys.stdin, help="Input file name")
#@click.option("-o", "--output", "outfile", type=click.File("w"), default=sys.stdout, help="Output file name")
#@click.option("-k", "--key", default="key", show_default=True, help="Sorting key")
@click.option('-n',        '--numberstart',         default = 2,          type=int,   help='Number of different starting genotypes.')		
@click.option('-a',        '--failure_a',           default = 24/(20/60), type=float, help='Failure rate of male gametocytes (per day).')		
@click.option('-b',        '--failure_b',           default = 24/(25/60), type=float, help='Failure rate of female gametocytes (per day).')			
@click.option('-r',        '--fertilization',       default = 0.08,       type=float, help='Fertilization of male and female gametes (per day).')	
@click.option('-i',        '--zygote_death',        default = 1,          type=int,   help='Death rate of zygotes (per day).')	
@click.option('-q',        '--transform_zygotes',   default = 24/19,      type=float, help='Transformation rate of zygotes (per day).')		
@click.option('-e',        '--transform_ookinetes', default = 0.6,        type=float, help='Transformation rate of ookinetes (per day).')		
@click.option('-j',        '--death_ookinetes',     default = 1.4,        type=float, help='Death rate of ookinetes (per day).')		
@click.option('-k',        '--death_oocysts',       default = 0,          type=float, help='Death rate of oocysts (per day).')		
@click.option('-o',        '--num_sporo',           default = 3e3,        type=int,   help='Number of sporozoites per oocyst.')		
@click.option('-p',        '--prop_sporo',          default = 0.2,        type=float, help='Proportion of sporozoites that make it to salivary gland.')		
@click.option('-w',        '--fraction_male',       default = 0.39,       type=float, help='Fraction of male gametes that are viable.')	
@click.option('-x',        '--fraction_female',     default = 1, 	      type=float, help='Fraction of female gametes that are viable.')	
@click.option('-y',        '--gametes_fem',         default = 1,          type=float, help='Number of female gametes per female gametocyte.')		
@click.option('-z',        '--gametes_male',        default = 8,	      type=int,   help='Number of male gametes per male gametocytes.')	
@click.option('-u',        '--k_param',             default = 1/7,        type=float, help='Parameter k.')	
@click.option('-m',        '--initial_time',        default = 10,         type=int,   help='Initial time value.')	
@click.option('-l',        '--max_bias',            default = 0.1,        type=float, help='The max bias value.')	
@click.option('-c',        '--number_sim',          default = 1,          type=int,   help='Number of simulations to run.')	
def main(numberstart, failure_a, failure_b, fertilization, zygote_death, transform_zygotes, transform_ookinetes, death_ookinetes, death_oocysts, num_sporo, prop_sporo, fraction_male, fraction_female, gametes_fem, gametes_male, k_param, initial_time, max_bias, number_sim):
    print("Running simulations(s)!")
    #Main routine code will go here. 

#param.N = 2; % Number of different starting genotypes
#param.a = 24/(20/60); % failure rate of male gametocytes (per day)
#param.b = 24/(25/60); % failure rate of female gametocytes (per day)
#param.r = 0.08; % fertilization of male and female gametes (per day)
#param.mu_z = 1; % death rate of zygotes (per day)
#param.sigma_z = 24/19; % transformation rate of zygotes (per day)
#param.sigma_e = 0.6; % transformation rate of ookinetes (per day)
#param.mu_e = 1.4; % death rate of ookinetes (per day)
#param.mu_o = 0; % death rate of oocysts (per day) **changed from paper**
#param.n = 3e3; % number of sporozoites per oocyst
#param.p = 0.2; % proportion of sporozoites that make it to salivary gland
#param.alpha = 0.39; % fraction of male gametes that are viable
#param.beta = 1;%0.96; % fraction of female gametes that are viable
#param.eta = 1; % number of female gametes per female gametocyte
#param.rho = 8; % number of male gametes per male gametocytes
#param.k = 1/7;
#param.t0 = 10;
#param.max_bias = 0.1;%0.5;%0.1;
#param.NumSim = 1;%10000;
    
#Example code concerning sub-commands:
    
#@cli.command()
#@click.argument('subcommand')
#@click.pass_context
#def help(ctx, subcommand):
#    subcommand_obj = cli.get_command(ctx, subcommand)
#    if subcommand_obj is None:
#        click.echo("I don't know that command.")
#    else:
#        click.echo(subcommand_obj.get_help(ctx))
        
if __name__ == "__main__":
    main()
    