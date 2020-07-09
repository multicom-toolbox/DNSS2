
###Evaluate DNSS2 for 3-class
```
perl lib/evaluation_ss_prediction.pl  -indir methods_prediction/DNSS2_results/ -list datasets/dnss2-blind-test.lst -ssa datasets/seq_labels -server "DNSS2" -tag DNSS2_test_DNSS2_prediction  -eva_ext '.ss_DNSS'  -seq_ext '.fasta'

#Total Matches: 78227   Total residues: 92858
#Total Q3 score    : 84.2436839044563
#Average Q3 score  : 84.6381116728131
#Average Sov score : 75.572052107119
```


###Evaluate DNSS2 for 8-class
```
perl lib/evaluation_ss_prediction_ss8.pl  -indir methods_prediction/DNSS2_results/ -list datasets/dnss2-blind-test.lst -ssa datasets/seq_labels -server "DNSS2" -tag DNSS2_test_DNSS2-ss8_prediction  -eva_ext '.ss8_DNSS'  -seq_ext '.fasta'

#Total Matches: 69165   Total residues: 92858
#Total Q8 score    : 74.4846970643348
#Average Q8 score  : 75.4631690193353
#Average Sov score : 75.504306636193
```

###Examples to evaluate DeepCNF, MUFOLD, psipred, PSSPred, SSPRO
```
perl lib/evaluation_ss_prediction.pl  -indir methods_prediction/DeepCNF_results/ -list datasets/dnss2-blind-test.lst -ssa datasets/seq_labels -server "DeepCNF" -tag DNSS2_test_DeepCNF_prediction  -eva_ext '.ss_DeepCNF'  -seq_ext '.fasta'

#Total Matches: 76215   Total residues: 92858
#Total Q3 score    : 82.0769346744491
#Average Q3 score  : 82.7555761192569
#Average Sov score : 70.2460603517276
```



###Model Training

###3-class training



###8-class training
