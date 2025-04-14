# script made for going through monophyly of a set of trees and outputting branch lengths
# made for core/peripheral phylogenetic simulations
# usage example: python summarize_branch.py -i 'sim_tree_output/*treefile' -o coresim_output_branch.csv

# ete3 main package used
import ete3, glob, argparse, re, pandas as pd

def main():
    # set up parser for output file
    parse = argparse.ArgumentParser(description = "Get output file")
    parse.add_argument("-o", "--output", type=str, help="Path to output file")
    parse.add_argument("-i", "--intrees", type=str, help="String for glob to interpret for trees, e.g. '/path/to/input/*tree' ")
    args = parse.parse_args()
    intrees, output = args.intrees, args.output
    
    # get all file names mathcing intrees string
    tree_names = glob.glob(intrees)
    meta_dict = {}
    # loops through trees
    for file in tree_names:
        tree = ete3.Tree(file)
        # split name, remove replicate number
        name = file.split('/')[2].split('_r')[0]
        # if name not in the dictionary, add it as a key
        if name not in meta_dict:
            # initialize with all variables, note that topo-specific distances not used in paper
            meta_dict[name] = {'Anest':0, 'Bnest':0, 'ABnest':0, 'Cmono':0, 'Other':0, 'Non-monophyletic':0, 'Poor-resolution':0, 'Anest_total':0, 'Bnest_total':0, 'ABnest_total':0, 'Cmono_total':0, 'Other_total':0, 'Non-monophyletic_total':0, 'Poor-resolution_total':0, 'Anest_core':0, 'Bnest_core':0, 'ABnest_core':0, 'Cmono_core':0, 'Other_core':0, 'Non-monophyletic_core':0, 'Poor-resolution_core':0, 'Anest_a':0, 'Bnest_a':0, 'ABnest_a':0, 'Cmono_a':0, 'Other_a':0, 'Non-monophyletic_a':0, 'Poor-resolution_a':0, 'Anest_b':0, 'Bnest_b':0, 'ABnest_b':0, 'Cmono_b':0, 'Other_b':0, 'Non-monophyletic_b':0, 'Poor-resolution_b':0, 'total':0, 'core':0, 'a':0, 'b':0}
            meta_dict[name]['name'] = name
            # parses parameters from name, notably removes periods, so decimals are taken out
            params = re.sub('\.', '', re.sub('[a-z]', '',  file)).split('_')
            # history, normally cut out of characters so include it separately
            hist = name.split('_')[5]
            meta_dict[name].update({'disp':params[1], 'ne':params[2], 'time':params[3], 'int':params[4], 'history':hist})
        # restrict to one outgroup sample and root
        tree.prune(['Outgroup1', 'A_1', 'A_2', 'A_3', 'A_4', 'B_1', 'B_2', 'B_3', 'B_4', 'Core1_1', 'Core1_2', 'Core1_3', 'Core1_4', 'Core2_1', 'Core2_2', 'Core2_3', 'Core2_4', 'Core3_1', 'Core3_2', 'Core3_3', 'Core3_4', 'Core4_1', 'Core4_2', 'Core4_3', 'Core4_4'], preserve_branch_length=True)
        tree.set_outgroup("Outgroup1")
        # only proceeds if A and B are monophyletic, doesn't care about core
        if check_mono(tree):
            tree.prune(['Outgroup1', 'A_1', 'B_1', 'Core1_1', 'Core2_1', 'Core3_1', 'Core4_1'], preserve_branch_length=True)
            # "get_topology" includes a lot of stuff, only want topology first
            topo_out = get_topology(tree)
            meta_dict[name][topo_out[0]]+=1 # increase count
            meta_dict[name][topo_out[0]+"_total"]+=topo_out[1] # add to topo-specific tree height sum
            meta_dict[name][topo_out[0]+"_core"]+=topo_out[2] # add to topo's mean core branch length sum
            meta_dict[name][topo_out[0]+"_a"]+=topo_out[3] # add to topo's a branch length sum
            meta_dict[name][topo_out[0]+"_b"]+=topo_out[4] # add to topo's b branch length sum
            meta_dict[name]["total"]+=topo_out[1] # add to parameter-combo-specific tree height sum
            meta_dict[name]["core"]+=topo_out[2] # add to parameter-combo-specific mean core branch length sum
            meta_dict[name]["a"]+=topo_out[3] # add to parameter-combo-specific a branch length sum
            meta_dict[name]["b"]+=topo_out[4] # add to parameter-combo-specific b branch length height sum

        else:
            meta_dict[name]['Non-monophyletic']+=1 # increase count for non-monophyletic
    # convert to data frame and output
    df = pd.DataFrame(meta_dict)
    df.transpose().to_csv(output)

# checks monophyley of A and B
def check_mono(tree):
    if tree.check_monophyly(['A_1','A_2','A_3','A_4'], target_attr="name")[0] & \
       tree.check_monophyly(['B_1','B_2','B_3','B_4'], target_attr="name")[0]:
        return True
    else:
        return False
        
def get_topology(tree):
    # subsets for stricting to first sample per population
    a = ['A_1']
    b = ['B_1']
    c = ['Core1_1', 'Core2_1', 'Core3_1', 'Core4_1']
    # makes a copy for assessing support among the 3 groups (just one core for now)
    supp_tree = tree.copy()
    supp_tree.prune(['Outgroup1', 'A_1', 'B_1', 'Core1_1'], preserve_branch_length=True)
    min_support=100.0
    # loops over support tree to play into output
    for node in supp_tree.traverse():
        # Support greater than 1.0 is to ignore leaves
        if(node.support<min_support and node.support>1.0): 
            min_support=node.support
    # loops over all ingroup to get the maximum and mean distance
    ingroup = a+b+c
    total_height = 0; n_t=0; max_height = 0;
    for i in ingroup:
        for j in ingroup:
            if i != j:
                dist = tree.get_distance(i, j)
                total_height  += dist
                n_t += 1
                if dist > max_height:
                    max_height = dist
    mean_height = total_height/n_t # double-counts but it gets averaged out
    # loops over distances between core samples
    total_core = 0
    for x in c:
        total_core += (tree&x).dist
    # sests mean core and a/b distances
    mean_core = total_core/4
    a_height = (tree&'A_1').dist
    b_height = (tree&'B_1').dist
    # if support is low, check is the core pops are well supported and monophyletic
    if min_support<70.0:
        if tree.check_monophyly(c, target_attr="name")[0]:
            c_tree = tree.copy()
            min_c = 100.0
            c_tree.prune(['Outgroup1', 'Core1_1', 'Core2_1', 'Core3_1', 'Core4_1'], preserve_branch_length=True)
            for node in c_tree.traverse():
                # Support greater than 1.0 is to ignore leaves
                if(node.support<min_c and node.support>1.0):
                    min_c=node.support
            # if core is monophyletic and well supported, note it as such regardless of other relationships
            if min_c >= 70.0:
                return(["Cmono",max_height, mean_core, a_height, b_height])
            # if not, mark as poorly resolved
            else:
                return(["Poor-resolution",max_height, mean_core, a_height, b_height])
        # if core is non-monophyletic and tree is otherwise poorly supported, make as poorly resolved
        else:
            return(["Poor-resolution",max_height, mean_core, a_height, b_height])
    # checks monophyly of certain groups to get at relevant topology
    elif tree.check_monophyly(c, target_attr="name")[0]:
        return(["Cmono",max_height, mean_core, a_height, b_height])
    elif tree.check_monophyly(a + c, target_attr="name")[0]:
        return(["Anest",max_height, mean_core, a_height, b_height])
    elif tree.check_monophyly(b + c, target_attr="name")[0]:
        return(["Bnest",max_height, mean_core, a_height, b_height])
    elif tree.check_monophyly(a + b + c, target_attr="name")[0]:
        return(["ABnest",max_height, mean_core, a_height, b_height])
    else:
        return(["Other",max_height, mean_core, a_height, b_height])
    # catch all, never hit in our simulations
    return(["Other",max_height, mean_core, a_height, b_height])


if __name__ == '__main__':
    main()
