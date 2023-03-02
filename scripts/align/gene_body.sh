mkdir -p Map/QC/geneBody_coverage
ls Map/*.bam |awk -F "." '{print $1}' | awk -F "/" '{print $2}' | while read id;do
	echo "/nfs1/public2/User/wzc/miniconda3/envs/rnaseq/bin/geneBody_coverage.py -r $(pwd)/Map/QC/ref.bed -i $(pwd)/Map/$id.bam -o $(pwd)/Map/QC/geneBody_coverage/$id.geneBodyCoverage.r";
done > geneBodyCoverage.sh
job array -s geneBodyCoverage.sh -n geneBody_$(date "+%Y%m%d%H%M%S") -m 50 -w
ls $(pwd)/Map/QC/geneBody_coverage/* | awk -F "geneBodyCoverage.r." '{print "mv " $0" "$1$2}' | bash
sed -i "s/geneBodyCoverage.r.//g" `ls $(pwd)/Map/QC/geneBody_coverage/*.r`
ls $(pwd)/Map/QC/geneBody_coverage/*.r | while read i;do
	echo "/nfs1/public2/User/LuRX/miniconda3/envs/rnaseq/bin/Rscript $i" | bash
done
ls $(pwd)/Map/QC/geneBody_coverage/*.pdf | awk -F ".pdf" '{print $1}' |while read i;do
	echo "/nfs1/public2/User/JXQ/miniconda3/bin/pdftoppm $i.pdf -singlefile -f 1 -png > $i.png" | bash
done

# /nfs1/public2/User/wzc/miniconda3/envs/rnaseq/bin/samtools index -@ 10 $(pwd)/Map/$id.bam $(pwd)/Map/$id.bam.bai && 


#  && mv $(pwd)/Map/QC/geneBody_coverage/$id.geneBodyCoverage.r.geneBodyCoverage.r $(pwd)/Map/QC/geneBody_coverage/$id.geneBodyCoverage.r
# $id 是sample名称