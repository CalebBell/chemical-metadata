cd data/crawl
xargs -d '\n' -a ../know_common_chemistry_CASRNs.txt -I arg -P4 -n1 bash -c 'wget "http://www.commonchemistry.org/ChemicalDetail.aspx?ref="arg'
