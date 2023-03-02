# -*- coding: utf-8 -*-
# 拆分步骤
"""
1.kegg_enrich文件结合diff_anno表格 提取每一条通路的全部KO
2.对每个KO，获得全部注释的gene以及上下调信息，全部基因上调为红色，全部基因下调为绿色，都有为蓝色
3.对于每个KO获取其在通路图中的位置信息
4.构建URL请求，并下载图片
5.爬取网页中的map部分的代码，更换注释信息
6.按照结构生成html，并保存
7.在表格中添加本地文件的超链接
"""
import re
import os
import sys
import random
import urllib
import requests
import pandas as pd
from bs4 import BeautifulSoup
from pathlib import Path


def splitKO(KOinfo: str):
    if not str(KOinfo) == "nan":
       return KOinfo.split(',')[0]
    return KOinfo


def getblockWeb(ko: str,KO_list: list)-> dict:
    """
    get the block html ,find KO - point(like 1.14.16.1 ...)        
    return dict of KO to point

    """
    KO_Point_dict = {}
    orignal_url = "https://www.kegg.jp/kegg-bin/show_pathway?"
    myid = re.sub(r'^ko','map',ko)
    # print(myid)
    blockurl = orignal_url+myid
    # print(blockurl)
    user_agent=['Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.87 Safari/537.36',
        'Mozilla/5.0 (X11; U; Linux x86_64; zh-CN; rv:1.9.2.10) Gecko/20100922 Ubuntu/10.10 (maverick) Firefox/3.6.10',
        'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/30.0.1599.101 Safari/537.36',
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.71 Safari/537.1 LBBROWSER',
        'Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E; QQBrowser/7.0.3698.400)',
        ]
    headers={
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Encoding': 'gzip, deflate, sdch',
        'Accept-Language': 'zh-CN,zh;q=0.8',
        'User-Agent': user_agent[random.randint(0,5)]
        }
    
    
    response = requests.get(blockurl,headers=headers)
    if response.text == '<html><body>'+myid+' not found.</body></html>':
        print("There is no {} ({})".format(ko,blockurl))
        return 'no kegg web'
    soup=BeautifulSoup(response.text,'lxml')
    
    # 解析全部的标签内容
    for content in soup.find(id ='mapdata').contents:
        content = str(content)

        if content is "\n" or content == "\r" or  content == "":
            next

        try:
            title = re.search(r".*title=\"(?P<title>.*)\"",content)
            title = title.group('title')


            ID = re.search(r".*\shref=\"(?P<ID>[\S]+)\"",content)
            ID = ID.group('ID')

            coords = re.search(r".*data-coords=\"(?P<coords>[\d,]+)\"",content)
            coords = coords.group('coords')
            
            refind_KOlist = re.findall(r'(K\d+)',title)
            refind_KOnumber = re.search(r',\s(?P<number>\d.\d.+), ',title)
            
            if len(refind_KOlist) == 0:
                next
            else:
                for reKO in refind_KOlist:
                    KO_Point_dict[reKO] = refind_KOnumber
        except:
            pass
            #print("Error om {} !".format(content))


    for content in soup.find(id ='module_mapdata').contents:
        content = str(content)

        if content is "\n" or content == "\r" or  content == "":
            next

        try:
            title = re.search(r".*title=\"(?P<title>.*)\"",content)
            title = title.group('title')


            ID = re.search(r".*\shref=\"(?P<ID>[\S]+)\"",content)
            ID = ID.group('ID')

            coords = re.search(r".*data-coords=\"(?P<coords>[\d,]+)\"",content)
            coords = coords.group('coords')
            
            refind_KOlist = re.findall(r'(K\d+)',title)
            refind_KOnumber = re.search(r',\s(?P<number>\d.\d.+), ',title)
            
            if len(refind_KOlist) == 0:
                next
            else:
                for reKO in refind_KOlist:
                    KO_Point_dict[reKO] = refind_KOnumber
        except:
            pass
            #print("Error om {} !".format(content))


        return KO_Point_dict


    
def getko_url(KO_Point_dict: dict,ko: str, KO_list: list, KO_dict: dict)-> str:
        """
        get the url to kegg database for colorful annotations of KOs and diff levels
        return a path of url
        e.g.
        https://www.kegg.jp/kegg-bin/show_pathway?map=map00100&multi_query=K07436+red,blue
        """
        #map00100&multi_query=K07436+red,blue
        part_url = 'https://www.kegg.jp/kegg-bin/show_pathway?map='
        ko_url = part_url + ko + '&multi_query='
        for KO in KO_list:
            level = set(KO_dict[KO]['level'].tolist()) 
#             print("level:")
#             print(level)
            level_color = 'blue,black'
            if level == set(['Increased']) or level == set(['up']) or level == set(['Up']):
                level_color = 'red,black'
            elif level == set(["Decreased"]) or level == set(['down']) or level == set(['Down']):
                level_color = 'green,black'
                
            ko_url = ko_url+KO+"+"+level_color+"%0d%0a"
        return ko_url
    
def getWeb(ko: str,KO_url: str,result_dir: str) -> list:
        """
        get the passway picture and passway html of kegg database.        
        return the ko_html,ko_pic  path
            
        """
        ko_html,ko_pic = 0,0
        work_dir = Path(result_dir)
        if not work_dir.exists():
            work_dir.mkdir(exist_ok=True)
            
        user_agent=['Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.87 Safari/537.36',
        'Mozilla/5.0 (X11; U; Linux x86_64; zh-CN; rv:1.9.2.10) Gecko/20100922 Ubuntu/10.10 (maverick) Firefox/3.6.10',
        'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/30.0.1599.101 Safari/537.36',
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.71 Safari/537.1 LBBROWSER',
        'Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E; QQBrowser/7.0.3698.400)',
        ]
        headers={'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Encoding': 'gzip, deflate, sdch',
        'Accept-Language': 'zh-CN,zh;q=0.8',
        'User-Agent': user_agent[random.randint(0,5)]
        }    
            
            
            
        response = requests.get(KO_url,headers=headers)
        
        # ------ 获取帖子内所有图片地址的方法 ------
        
        soup=BeautifulSoup(response.text,'lxml')
        imginfo = soup.find('img',id="pathwayimage")
        imgsrc = str(re.search(r"src=\"(?P<src>.*\S*\.png)\"", str(imginfo)))

        imgsrc_path = Path(Path(imgsrc).name).stem
        imgName = imgsrc_path

        # ------ 这里最好使用异常处理及多线程编程方式 ------
        pattern = r'src="(.*\S*\.png)"'
        # try:      
        imgsrc_path = re.search(pattern, str(imginfo)).group(1)
        # catch AttributeError：
            
        try:


            img_data = urllib.request.urlopen('https://www.kegg.jp'+imgsrc_path).read()
            ko_pic = os.path.join(result_dir , str(imgName)+".png")

            f = open(ko_pic, 'wb')
            f.write(img_data)
            f.close()
        except Exception as e:
            print('https://www.kegg.jp/'+imgsrc_path+ "error")
            
        htl = open(os.path.join(result_dir , str(imgName)+"_old.html"),'w')
        htl.write(soup.prettify())
        htl.close()
        ko_html = os.path.join(result_dir ,str(imgName)+"_old.html")

        return ko_html,ko_pic
"""
# 测试重构web
1.Beautiful 提取全部area ，并转化，过滤  
    1.提取样式信息  
    2.将title结合KO_dict[KOID]，结合level以及logFC信息，转化成info  
2.重建web code  
"""
def creatNew(ko:str ,ko_html: str,ko_pic: str,KO_list: list,KO_dict: dict)-> str:
        """
        read the html , build a new html  
        return the path of new html
            
        """
        new_dir = os.path.dirname(ko_html)
        path_html = os.path.join(new_dir,ko+'_new.html')
        print(ko,ko_html,ko_pic,KO_list)
        
        # 解析ko_html
        ko_html_file = open(ko_html,'r')
        ko_html_text = ko_html_file.read()
        #print(ko_html_text)
        soup = BeautifulSoup(ko_html_text)
        #print(soup)
        area_list = soup.find_all('area')
        
        # html 框架
        html_head = "<html>\n<head>\n<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n<title>"
        html_head = html_head + ko
        html_head = html_head +"""</title>
<style type="text/css">
#result ul{
    width:95%;
    list-style: none;
}
#result ul li{
    display:inline;
    word-break:break-all;
    word-wrap : break-word ;
}
</style>
<script type="text/javascript">
function showInfo(info) {
    obj = document.getElementById("result");
    obj.innerHTML = "<div style='cursor: pointer; position: absolute; right: 5px; color: #000;' onclick='javascript: document.getElementById(\\\"result\\\").style.display = \\\"none\\\";' title='close'>X</div>\" + info;
    obj.style.top = document.body.scrollTop;
    obj.style.left = document.body.scrollLeft;
    obj.style.display = "";
}
</script>
</head>
<body>
<map name="""+ko+">\n"
        html_tail = """</map>
<img src='"""+ko_pic+"""' usemap='#"""+ko+"""' />
<div id='result' style='position: absolute; width: 50%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;' onmouseover="javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;" onmouseout="javascript: this.style.filter = 'alpha(opacity=95)'; this.style.opacity = 0.95;"></div>
</body></html>"""
    
        
        #整理全部
        new_html = open(path_html,'w') 
        new_html.write(html_head)
        
        for area in area_list:
#             print('Before...')
#             print(area)
            try:
                area['href'] = "http://www.kegg.jp/"+ area['href']
            except:
                error_info = "Cant create {}".format(path_html)
                continue

            try:
                area['onmouseover'] = "\n"  + Title2JsInfo(area['title'],KO_dict)
                del area['data-coords']
                del area['data-entry']
                del area['title']
                del area['class']
            except:
                area['onmouseover'] = 'None'
                del area['data-coords']
                del area['data-entry']
                del area['class']
#             print('After...')
#             print(area)
            if area['onmouseover'] == "\n"+"None informations" or area['onmouseover'] == 'None':
                del area['onmouseover']
                continue
            new_html.write(area.prettify()+'\n')
            
#             print(area.prettify())
            #print(area.encode("utf-8"))
        new_html.write(html_tail)
        new_html.close()
        return path_html

def Title2JsInfo(title_info:str,KO_dict: dict)-> str:
    """
    get info about diff genes from title 
    return to the code of js ,to show some infos (Up/Down,geneID,logFC..)
    
    """
#     print(title_info)
    ID_list = title_info.split(',')
    ID_list2 = []
    for temp_id in ID_list:
        temp_id = temp_id.strip()
#        print(temp_id)
        ID_list2.append(temp_id.split(' ')[0])
#     print(ID_list2)
    
    Up_info = ''
    Down_info = ''
    Up_KO = []
    Down_KO = []
        
    for ID in ID_list2:
        
        if ID in KO_dict.keys():
            myKO_dict = KO_dict[ID]
            UP_df = myKO_dict[myKO_dict['level'] == 'Increased'] #or myKO_dict['level'] == 'up' or myKO_dict['level'] == 'UP']]
            Down_df = myKO_dict[myKO_dict['level'] == 'Decreased']# or myKO_dict['level'] == 'down' or myKO_dict['level'] == 'Down' ]

            if len(UP_df) >= 1:
                Up_info = KODictDF2str(UP_df)
                Up_KO_list = list(set(UP_df['KOID'].tolist()))
                if Up_KO == []:
                    Up_KO = Up_KO_list
                else:
                    Up_KO = Up_KO +Up_KO_list
            if len(Down_df) >= 1:
                Down_info = KODictDF2str(Down_df)
                Down_KO_list = list(set(Down_df['KOID'].tolist()))
                if Down_KO == []:
                    Down_KO = Down_KO_list
                else:
                    Down_KO = Down_KO +Down_KO_list
                
    if Up_KO != [] and Down_KO != []:
        jscode = "javascript: showInfo(\"<ul><li style=\\\"color: #f00;\\\">"+ title_info + ':'+"Up regulated<ul><li>"+Up_info+"</li></ul></li></ul><ul><li style=\\\"color: #0f0;\\\">Down regulated<ul><li>"+Down_info+"</li></ul></li></ul>\");"#.encode("utf-8")
    elif Up_KO == [] and Down_KO == []:
        jscode = "None informations"
    elif Up_KO != [] and Down_KO == []:
        jscode = "javascript: showInfo(\"<ul><li style=\\\"color: #f00;\\\">"+ title_info + ':'+"Up regulated<ul><li>"+Up_info+"</li></ul></li></ul>\");"#.encode("utf-8")
    elif Up_KO == [] and Down_KO != []:
        jscode = "javascript: showInfo(\"<ul><li style=\\\"color: #0f0;\\\">"+ title_info + ':'+"Down regulated<ul><li>"+Down_info+"</li></ul></li></ul>\");"#.encode("utf-8")

    return jscode

def KODictDF2str(df)-> str:
    '''df:
    gene_id    KOID    level    logFC
670    PH01000295G0710    K00006    Increased    2.368777
3547    PH01001609G0120    K00006    Decreased    -1.593535
4282    PH01000393G0300    K00006    Decreased    -2.086649
    '''
    info_all = ''
    for indexs in df.index:
        info = df.loc[indexs,['gene_id','logFC']].values
#         print(info.tolist())

        info = str(info[0])+'('+str(info[1])+')'
        if info_all == '':
            info_all = info
        else:
            info_all = info_all+','+info
    return info_all

def main(diff_anno_file,kegg_anno,result_dir):
    #1.kegg_enrich文件结合diff_anno表格 提取每一条通路的全部KO
    # diff_anno = pd.read_csv('B2_vs_B1_diff_anno.xls',sep="\t") 
    # KO_Group = gene2KO.groupby("KOID")

    diff_anno = pd.read_csv(diff_anno_file,sep="\t") 
    

    diff_anno['KOID'] = diff_anno['KO'].apply(splitKO)
    diff_anno = diff_anno.dropna(axis=0)

    gene2KO = diff_anno.loc[:,['gene_id','KOID','level','logFC']]
    gene2KO = gene2KO[gene2KO['level'] != 'nonsignificant']

    KO_Group = gene2KO.groupby("KOID")
    # 2.对每个KO，获得全部注释的gene以及上下调信息，全部基因上调为红色，全部基因下调为绿色，都有为蓝色
    KO_dict = {}
    for OneKO in KO_Group:
        KO_dict[OneKO[0]] = OneKO[1]

    # enrich_df = pd.read_csv('kegg_enrich_all.xls',sep="\t") 
    enrich_df = pd.read_csv(kegg_anno,sep="\t")

    for ko in enrich_df['ID'].tolist():
        genes = enrich_df[enrich_df['ID'] == ko]['geneID'].tolist()
        try:
            genestr = genes[0]
        except IndexError:
            genes = str(genes)
            if "'" in genes:
                genestr = genes.split("'")
            else:
                genestr = genes.split("\"")
        print(genestr)
        gene_list = genestr.split('/')
        gene_list = list(set(gene_list))
        KO_list = []
        for gene in gene_list:
            mylist = gene2KO[gene2KO['gene_id'] == gene]["KOID"].tolist()
            if(len(mylist)!=0):
                for i in mylist:
                    KO_list.append(i)
        KO_list = list(set(KO_list))
        # print(KO_list)
        # 1，获得无颜色界面，然后解析html页面，根据标注找到关键KO的点（类似1.14.16.1）
        KO_Point_dict = getblockWeb(ko, KO_list)
        if KO_Point_dict == 'no kegg web':
            continue
        # 2，根据KO-点-上下调，获得颜色注释界面，并下载html以及图片
        ko_url = getko_url(KO_Point_dict, ko, KO_list, KO_dict)
        # html_path,img_path = getWeb(ko, ko_url, 'test')
        try:
            html_path,img_path = getWeb(ko, ko_url, result_dir)
        except AttributeError:
            continue
        print(html_path,img_path)
        # 3，根据上下调，重构html，添加注释效果
        # 获得图片相对于页面的相对路径
        try:
            img_dir ,img_name = os.path.split(img_path)
            img_relative_path = './'+img_name
            getnew_html = creatNew(ko,html_path, img_relative_path, KO_list, KO_dict)
        except TypeError:
            print("Error HTML: " + ko_url)
        # print(getnew_html)


if __name__ == '__main__':
    args_list = sys.argv[1:]

    diff_anno = args_list[0]
    kegg_anno = args_list[1]
    result_dir = args_list[2]

    main(diff_anno,kegg_anno,result_dir)

