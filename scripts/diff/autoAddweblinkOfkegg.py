#conding:utf-8
import pandas as pd
import requests
import urllib
import os,sys,re
import argparse
# ======sub function =============
#
def keggAPI(ko_ID,ko_dict,ipath_result_dir):
	url = 'https://pathways.embl.de/mapping.cgi'
	headers = {'Content-Type':"application/x-www-form-urlencoded",
			   'Cookie':"_pk_ses.7.a35f=1; _ga=GA1.2.883760841.1593483374; _gid=GA1.2.550971889.1593483374; _pk_id.7.a35f=8f615d9cdfff502f.1593482137.1.1593485705.1593482137.",
			  'User-Agent':"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.116 Safari/537.36",
			   'Accept':"svg/html,text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9"}
 
	KO_list = ko_dict[ko_ID]
	KO_str = '%0D'.join(KO_list) 

	data = "selection=" + KO_str + "&default_opacity=1&default_width=3&default_radius=7&default_color=%23aaaaaa&background_color=%23ffffff&tax_filter=&map=metabolic&export_type=svg&export_dpi=120"
	outfile = os.path.join(ipath_result_dir,ko_ID+".html")
	r = requests.post(url,data,headers)
	r_xml = r.text
	f = open(outfile,'w') 
	f.write(r_xml)
	f.close()

def g2K(g,K):
	gene_list = g.to_list()
	Ko_list = K.to_list()
	res_dic = {}
	
	for i in range(len(gene_list)):
		res_dic[gene_list[i]] = Ko_list[i]
		
	return res_dic

def k2gs(ko,genes):
	ko_list = ko.to_list()
	genes_list = genes.to_list()
	res_dic = {}
	
	for i in range(len(ko_list)):
		res_dic[ko_list[i]] = genes_list[i]	
	return res_dic

def changeDic(k2gs_dic,g2k_dic):
	res_dic = {}
	for k,v in k2gs_dic.items():
		
		genes = v.split('/')
		#sys.exit(genes)
		KO_list = []
		for i in genes:
			KO = g2k_dic[i]
			KO_list.append(KO)
		new_v = '%0D'.join(KO_list)
		res_dic[k] = new_v
	return res_dic

# 添加KEGG weblink
def getweblink(ko_id):
	weblink = "=HYPERLINK(\"https://www.kegg.jp/kegg-bin/show_pathway?{}\",\"{} weblink\") ".format(ko_id,ko_id)
	return weblink

# 添加本地的ipath3 link
def getlocallink(ko_id,res_dir):
	filepath = os.path.join(res_dir,ko_id+'_new.html')
	locallink = "=HYPERLINK(\"{}\",\"{} local link\")".format(filepath,ko_id)
	return locallink

def WeblinkPath(All,Out):
	if Out:
		all_dic,all_file = os.path.split(All)
		all_name,all_ext = os.path.splitext(all_file)
		all_newpath = os.path.join(Out,all_name+'_weblink'+all_ext)


		return all_newpath
	else:
		all_file,all_ext = os.path.splitext(All)
		all_newpath = all_file+'_weblink'+all_ext

		return all_newpath


#
# ================================

# ====== main ====================
#
def main(Anno2Gene,KeggEnrichAll,DiffAnno,Pics,OutputDic):

	kegg_anno_df = pd.read_csv(Anno2Gene,sep='\t')
	diff_df = pd.read_csv(DiffAnno,sep='\t')
	kegg_enrich_df = pd.read_csv(KeggEnrichAll,sep="\t")

	g2k_dic = g2K(kegg_anno_df["gene_id"],kegg_anno_df["ko"]) 
	k2gs_dic = k2gs(kegg_enrich_df['ID'],kegg_enrich_df['geneID'])
	ko2KO_dict = changeDic(k2gs_dic,g2k_dic)


	kegg_all_df = pd.read_csv(KeggEnrichAll,sep="\t")
	kegg_all_df['weblink'] = kegg_all_df["ID"].apply(getweblink)


	if not Pics:
		pic_dic = '../maps'
	else:
		pic_dic = Pics

	kegg_all_df['locallink'] = kegg_all_df["ID"].apply(getlocallink,args=[pic_dic])


	allOutFile = WeblinkPath(KeggEnrichAll,OutputDic)

	kegg_all_df.to_excel(allOutFile,index=False)

	# for key,value in ko2KO_dict.items():
	# 	keggAPI(key,ko2KO_dict,"test")


if __name__ == '__main__':
	# args_lst = sys.args[1:]

	# if len(args_lst) < 4:
	# 	help_info = 'help'
	# 	print(help_info)
	
	parser = argparse.ArgumentParser() 


	parser = argparse.ArgumentParser();
	parser.add_argument("-a", "--Anno2Gene",help="a kegg_anno_gene2 file")
	parser.add_argument("-k", "--KeggEnrichAll", help="diff/kegg_enrich/kegg_enrich_all.xls")
	
	parser.add_argument("-kd", "--KeggEnrichDown", help="diff/kegg_enrich/kegg_enrich_down.xls")
	parser.add_argument("-ku", "--KeggEnrichUp", help="diff/kegg_enrich/kegg_enrich_up.xls")

	parser.add_argument("-d", "--DiffAnno",help="diff/diff_diff_anno.xls")
	parser.add_argument("-p", "--Pics",help="kegg 通路图的存放路径 注意：这里要填相对于输出文件的路径")
	parser.add_argument("-o", "--OutputDic", action="store_true",default=False,help="默认输出到kegg富集目录下，若有其他指定目录则传入路径即可")
	args = parser.parse_args()


	help = """optional arguments:
-h, --help            show this help message and exit
-a, --Anno2Gene       a kegg_anno_gene2 file
-k, --KeggEnrichAll   diff/kegg_enrich/kegg_enrich_all.xls
-kd, --KeggEnrichDown diff/kegg_enrich/kegg_enrich_down.xls
-ku, --KeggEnrichUp   diff/kegg_enrich/kegg_enrich_up.xls
-d, --DiffAnno        M2_vs_M1/M2_vs_M1_diff_anno.xls
-p, --Pics            kegg 通路图的存放路径
-o, --OutputDic       默认输出到kegg富集目录下，若有其他指定目录则传入路径即可
	
	"""

	main(args.Anno2Gene,args.KeggEnrichAll,args.DiffAnno,args.Pics,args.OutputDic)

#
# ================================