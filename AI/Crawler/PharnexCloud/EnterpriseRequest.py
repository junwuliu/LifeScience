#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
import argparse
import random
import time
import sys
import re
import os
import json
import requests
import pandas as pd
bindir = os.path.abspath(os.path.dirname(__file__))
global session
session = requests.session()

import atexit 
@atexit.register 
def clean(): 
    loginout_url = "https://db-server.pharnexcloud.com/external/authenticate/authorization"
    session.delete(loginout_url)
    print('结束登录')

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-g','--gene',help='gene symbol',dest='gene',type=str,required=True)
	parser.add_argument('-o','--outdir',help='outdir file',dest='outdir',type=str,required=True)
	args=parser.parse_args()
	#session = requests.session()
	genelist =  pd.read_csv(args.gene,header=0)['gene']
	## 模拟浏览器
	session.headers = { 
        "User-Agent": "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/106.0.0.0 Safari/537.36 "
    }
	## 企业版登录网址
	login_url = "*****"
	login_data = { 
    'username': '***',
    'password': '***',
    'NVCVal' : "***"
	}
	rep = session.post(login_url, data=login_data) ## 执行登录
	#print (login_url,login_data)
	if (rep.status_code == 200):
		print("Login Sucess")
		time.sleep(random.randint(1,3))
	else:
		print("Login Failed:" + str(rep.status_code))
	response = json.loads(rep.text)
	auth = response['token_type'] + " " + response['access_token']
	session.headers.update({"authorization": auth }) ## 更新token到头文件
	for gene in genelist:
		Scrapy(session,gene,args.outdir)
		time.sleep(random.randint(2,5))
	loginout_url = "***login_URL***"
	session.delete(loginout_url)
	print ("Done")
	
def Scrapy(session,gene,outdir):
	#print (session.headers)
	print (gene)
	## 获取搜索的基本信息
	basic_url = "***basic_url/medical/search?q=" + gene + "&is_manual_input=1"
	basic_info = session.get(basic_url)
	#print (basic_info.text)
	basic_dict = json.loads(basic_info.text)
	basic_data = basic_dict['data']
	if (basic_data):
		basic = pd.DataFrame(basic_data[0]['databases'])
	else:
		print (gene, "have no result !")
		return
	## 遍历大类
	for num in range(0,len(basic_data)):
		catagory = basic_data[num]['category']
		#print (pd.DataFrame(basic_data[num]['databases']))
		if (num > 0):
			basic = pd.concat([basic,pd.DataFrame(basic_data[num]['databases'])],axis=0,join='outer').reset_index(drop=True)
		for dbnum in range(0,len(basic_data[num]['databases'])): ## 遍历小类
			time.sleep(random.randint(2,4))
			response_result = pd.DataFrame()
			db_id = basic_data[num]['databases'][dbnum]['database_id']
			db_name = basic_data[num]['databases'][dbnum]['name']
			db_total = basic_data[num]['databases'][dbnum]['total']
			max_page = int((db_total+0.1)/50) + 1 if db_total % 50 == 0 else int((db_total+0.1)/50) + 2
			print ("现在爬取:", db_name)
			if (db_id == 441 ):
				## 全球药物研发数据
				#print ("start to get drug info")
				#request_url = "https://db-server.pharnexcloud.com/external/customization/ca/441/aggregations/drug"
				request_url = "***drug_url"
				#id = 'results'
				id = 'data'
			else:
				request_url = "***url/"+str(db_id)+"/documents"
				id = 'data'
			for pagenum in range(1,max_page):
				#time.sleep(1)
				post_data = {  'query' : gene,'per_page':50,'page' :  pagenum}
				#print (post_data)
				info = session.post(request_url,data = post_data)
				if (info.status_code == 200): ## 看是否有访问权限
					tmp = json.loads(info.text)
					response_result = pd.DataFrame(tmp[id]) if pagenum == 1 else pd.concat([response_result,pd.DataFrame(tmp[id])],axis=0,join='outer').reset_index(drop=True)
				else:
					continue
			if (response_result.empty):
				pass
			else:
				file = outdir + "/" + gene + "." + db_name + ".xls"
			#response_result.to_csv(file,sep="\t",encoding='utf-8_sig')
				response_result.to_csv(file,sep="\t",encoding='gb18030')
	print (basic)
	basicfile = outdir +  "/" + gene + ".BasicSearch.tsv"
	basic.to_csv(basicfile,sep="\t",encoding='utf-8_sig')
	#loginout_url = "loginout_url"
	#session.delete(loginout_url) ## login out

if __name__ == '__main__':
	main()
