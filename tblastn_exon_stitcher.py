import re; import os; import Bio; from Bio import SeqIO
cwd=os.getcwd()
dir_w_alignments=cwd+'/test_tblastn_data/'
outfilename=cwd+'/test.faa'

with open(outfilename, 'a') as outfile:
    for aln_filename in os.listdir(dir_w_alignments):
        hypothetical_protein=''; exon= ''; new_seqid=''; list_of_used_midpoints=[]; list_of_exons_in_this_hit=[]; total_exons=0; query_length=0; seqcount = 0; start_list=[]; pos=0; exon_dict={}
        if aln_filename.endswith(".aln"):
        
            print('\n********************************\n\nProcessing ' + aln_filename)
            for record in SeqIO.parse(dir_w_alignments + aln_filename, "fasta"):
                seqcount +=1; exon=''; hypothetical_protein=''
                if seqcount == 1: query_length = len(record.seq)## If it's the query sequence, check it's length
                else:
                    list_of_exons_in_this_hit = re.sub("---+","_", str(record.seq)).split('_') ##I consider any AAs separated by more than two gaps as separate exons

                    ### calculate the midpoint position for each exon
                    for exon in list_of_exons_in_this_hit:
                        if exon != '': ###                                #if the exon is not empty, use it
                            start_pos=record.seq.find(exon)#              #find starting_position for exon
                            exon = exon.replace('*', '').replace('-','')  #remove stop codons and gaps
                            half_length=len(exon)/2                       #; print(half_length)
                            midpoint=start_pos+half_length                #; print(str(midpoint))
                            if midpoint not in list_of_used_midpoints:
                                exon_dict[midpoint]=exon                      #Save to dictionary {midpoint:exon}
                                list_of_used_midpoints.append(midpoint)
                            else:
                                print('This exon is redundant and wont be concatenated to any other exons')
                                new_seq_id=aln_filename.replace('.aln','')+'b_1exons'; hypothetical_protein=exon
                                outfile.write('>'+new_seq_id+'\n'+hypothetical_protein+'\n')
                                hypothetical_protein=''

            ### sort the exons in order by their midpoints
            exon_midpoints_ordered_list=sorted(exon_dict.keys())
            for key in exon_midpoints_ordered_list:
                hypothetical_protein = hypothetical_protein + exon_dict[key] ## concatenate exons in the right order
    
            ### How many exons in total?
            total_exons = len(exon_dict); print('Total exons:' + str(total_exons))
    
            ### check how long the hypothetical protein is compared to the query protein
            print('Query length: ' + str(query_length) + '\n' + 'Aligned AA length:' + str(len(hypothetical_protein)))

            ### print that seqid and protein
            new_seq_id=aln_filename.replace('.aln','')+'_'+str(total_exons)+'exons'
        
            outfile.write('>'+new_seq_id+'\n'+hypothetical_protein+'\n')

