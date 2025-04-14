'''
Made by Ethan Gyllenhaal (egyllenh@ttu.edu)
Last updated 14 April 2025

Script for running a single msprime simulation for four core and two peripheral populations
Takes in parameters for output path, replicate number, number of migrants per generation (generally <1), per population effective pop. size, time to run after colonization completes, interval of colonization events, and a string split order.
Outputs a fasta for a subset of 2 individuals (4 chromosomes) per population.

Useage like python $dir/core_trial_noanc.py -r [int replicate] -m [float Nm] -n [int/float effective pop size] -e [int/float post-split evolution time] -i [int/float split interval in gens]  -s [string split interval: CAB, ACB, or ABC] -o [output stem, ideally includes all variables used]

run_msprime is the method for running the simulation and calculating summary statistics
main is simply the driver for everything
'''

import os, sys, math, msprime as msp, numpy as np, re, argparse, tskit, allel

def main():
    
    # set up parser and arguments with ArgParse
    parse = argparse.ArgumentParser(description = "Get simulation parameters")
    parse.add_argument("-o", "--outstem", type=str, help="Path to output directory and stem")
    parse.add_argument("-r", "--rep", type=int, help="Simulation replicate")
    parse.add_argument("-m", "--migration", type=float, help="Number of migrants per generation among core population")
    parse.add_argument("-n", "--ne", type=float, help="Per-population pop size")
    parse.add_argument("-e", "--evo", type=float, help="Time to run after colonization is done")
    parse.add_argument("-i", "--interval", type=float, help="Interval of colonziation events")
    parse.add_argument("-s", "--splitorder", type=str, help="String of split order, expects CAB, ACB, or ABC")
    args = parse.parse_args()
   
    # assign arguments to variables
    outstem, repnum, mig, popsize, evo, interval, split  = args.outstem, args.rep, args.migration, args.ne, args.evo, args.interval, args.splitorder
    
    # calls msprime, assigns results to variables, and writes those to output
    muts = run_msprime(mig, popsize, evo, interval, split)
    muts.write_fasta(outstem + ".fasta")


def run_msprime(mig, size, evo, interval, split):
    # set number of samples (wind up with 2 fasta/sample)
    samples=2

    # Set up final populations, which are set to be inially active if something splits from them
    demography = msp.Demography()
    demography.add_population(name="C1", initial_size=size, initially_active=True)
    demography.add_population(name="C2", initial_size=size, initially_active=True)
    demography.add_population(name="C3", initial_size=size, initially_active=True)
    demography.add_population(name="C4", initial_size=size, initially_active=True)
    demography.add_population(name="A", initial_size=size)
    demography.add_population(name="B", initial_size=size)
    demography.add_population(name="Out", initial_size=10000, initially_active=True)
    
    # note that interval for the core is split into sub-intervals for splits
    
    # core -> peripheral a -> peripheral b
    if(split=="CAB"):
        # Split order
        demography.add_population_split(time=evo+(4*interval), derived=["C1"], ancestral="Out")
        demography.add_population_split(time=evo+(3*interval), derived=["C2"], ancestral="C1")
        demography.add_population_split(time=evo+(2*interval)+((2*interval)//3), derived=["C3"], ancestral="C2")
        demography.add_population_split(time=evo+(2*interval)+(interval//3), derived=["C4"], ancestral="C3")
        demography.add_population_split(time=evo+(2*interval), derived=["A"], ancestral="C2")
        demography.add_population_split(time=evo+interval, derived=["B"], ancestral="C4")
    # peripheral a -> core -> peripheral b    
    elif(split=="ACB"):
        # Split order
        demography.add_population_split(time=evo+(4*interval), derived=["C1"], ancestral="Out")
        demography.add_population_split(time=evo+(3*interval), derived=["A"], ancestral="C1")
        demography.add_population_split(time=evo+(2*interval), derived=["C2"], ancestral="C1")
        demography.add_population_split(time=evo+(interval)+((2*interval)//3), derived=["C3"], ancestral="C2")
        demography.add_population_split(time=evo+(interval)+(interval//3), derived=["C4"], ancestral="C3")
        demography.add_population_split(time=evo+interval, derived=["B"], ancestral="C4")
    # peripheral a -> peripheral b -> core   
    elif(split=="ABC"):

        # Split order
        demography.add_population_split(time=evo+(4*interval), derived=["C1"], ancestral="Out")
        demography.add_population_split(time=evo+(3*interval), derived=["A"], ancestral="C1")
        demography.add_population_split(time=evo+(2*interval), derived=["B"], ancestral="C1")
        demography.add_population_split(time=evo+(interval), derived=["C2"], ancestral="C1")
        demography.add_population_split(time=evo+((2*interval)//3), derived=["C3"], ancestral="C2")
        demography.add_population_split(time=evo+(interval//3), derived=["C4"], ancestral="C3")

    # No dynamism! So just set migration rates
    # Note that migration is in #/gen, so needs to be divided by forward sink (backwards source) pop size
    # Fortunately here size is constant
    demography.set_symmetric_migration_rate(["C1", "C2"], rate=mig/size)
    demography.set_symmetric_migration_rate(["C2", "C3"], rate=mig/size)
    demography.set_symmetric_migration_rate(["C3", "C4"], rate=mig/size)
    
    # sort the demographic events
    demography.sort_events()

    # Run the simulation to get tree sequences
    trees=msp.sim_ancestry(samples={"C1":samples, "C2":samples, "C3":samples, "C4":samples, "A":samples, "B":samples, "Out":samples}, 
                           demography=demography, 
                           recombination_rate=1e-7,
                           sequence_length=1e7)
                           
    # Add mutations to treeseqs and output
    mutation=msp.sim_mutations(trees, rate=4.6e-9, model="jc69")
    return(mutation)

if __name__ == '__main__':
    main()
    
