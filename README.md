Ce répertoire vise à étudier l'influence de la méthylation sur le binding des ARFs.

Les données ont étés téléchargées depuis :

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1444040

Les trois fichiers :

GSM1444040_Col0_NCBI_depth4_mCG.wig.gz  
GSM1444040_Col0_NCBI_depth4_mCHG.wig.gz  
GSM1444040_Col0_NCBI_depth4_mCHH.wig.gz  

ont étés décompressés, puis converties au format bedGraph  grace au programme "wigToBedGraph" disponible sur le site ucsc (du projet encode).


Les fichiers ont été respectivement renomés en :

mCG.bedGraph  
mCHG.bedGraph  
mCHH.bedGraph  


Pour comprendre le fonctionnement des fichiers, voici la réponse de Kai Tang qui a généré les données.

------------------------------------------------------------------

Negative value in wig file is for Cs in "-"(minus) strand.

For example, in "GSM1444040_Col0_NCBI_depth4_mCHG.wig",   
variableStep    chrom=chr1  
32    0.4  
34    -0.5  

 
this means  for chr1, the 32nd base is "C" and its methylation level is 40%. While 34th base is "G" and the methylation level for its complement "C" in minus strand is 50%.

------------------------------------------------------------------