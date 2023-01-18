
#macse
parallel --jobs 13 java -jar /home/kpotoh/tools/macse_v2.06.jar -prog alignSequences -gc_def 2 -seq {} -out_AA mulal/{/.}.faa -out_NT mulal/{/.}.fna ::: genes/*

#iqtree
parallel --jobs 13 iqtree2 -s {} -m GTR+FO+R6+I -nt 4 --prefix {/.} --rate ::: mulal/*.fna
