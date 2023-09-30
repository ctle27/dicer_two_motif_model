# dicer_two_motif_model
1.make_ref.py is used to make all substrate sequences to be used for mapping and determine cleavage sites from the sequencing of cleavage fragments

2.Calculate_efficiency_accuracy_score.py is used to process mapping data of control shRNA samples and cleavage samples (double cleavage and single cleavage) to get the RPM value for each variant/cleavage fragment. From that, the cleavage efficiency and accuracy for each cleavage site of each shRNA variants were calculated.

3.Motif_identification_and_plotting.py utilizes the data generated from 4.Calculate_efficiency_accuracy_score.py to identify the structural features and sequence motifs that show high cleavage accuracy/efficiency.
