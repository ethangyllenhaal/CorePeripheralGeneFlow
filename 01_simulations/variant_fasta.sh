#!/bin/bash
# simple script for processing and rename fastas, with output being variant only and non-interleaved
# Useage: sh variant_fasta.sh in.fasta out.fasta


# processes input, removing Ns, finding relevant lines, removing newlines for all but name
tr -d 'N' < $1 |
	grep '[A,G,T,C,>]' |
	tr -d '\n' |
	sed 's/>n[0-9][0-9]*/&\n/g' | 
	sed 's/>/\n>/g' | tail -n +2 > $2

# a bunch of seds to select relevant lines and rename
# order based on msprime output
sed -i 's/>n24/>Outgroup1/g' $2
sed -i 's/>n25/>Outgroup2/g' $2
sed -i 's/>n26/>Outgroup3/g' $2
sed -i 's/>n27/>Outgroup4/g' $2

sed -i 's/>n16/>A_1/g' $2
sed -i 's/>n17/>A_2/g' $2
sed -i 's/>n18/>A_3/g' $2
sed -i 's/>n19/>A_4/g' $2

sed -i 's/>n20/>B_1/g' $2
sed -i 's/>n21/>B_2/g' $2
sed -i 's/>n22/>B_3/g' $2
sed -i 's/>n23/>B_4/g' $2

sed -i 's/>n4/>Core2_1/g' $2
sed -i 's/>n5/>Core2_2/g' $2
sed -i 's/>n6/>Core2_3/g' $2
sed -i 's/>n7/>Core2_4/g' $2

sed -i 's/>n8/>Core3_1/g' $2
sed -i 's/>n9/>Core3_2/g' $2
sed -i 's/>n10/>Core3_3/g' $2
sed -i 's/>n11/>Core3_4/g' $2

sed -i 's/>n12/>Core4_1/g' $2
sed -i 's/>n13/>Core4_2/g' $2
sed -i 's/>n14/>Core4_3/g' $2
sed -i 's/>n15/>Core4_4/g' $2

sed -i 's/>n0/>Core1_1/g' $2
sed -i 's/>n1/>Core1_2/g' $2
sed -i 's/>n2/>Core1_3/g' $2
sed -i 's/>n3/>Core1_4/g' $2
