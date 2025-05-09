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
    loginout_url = "https://db-server/external/authenticate/authorization"
    session.delete(loginout_url)
    print('结束登录')

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-p','--patent',help='patent list',dest='patent',type=str,required=True)
	parser.add_argument('-o','--outdir',help='outdir file',dest='outdir',type=str,required=True)
	args=parser.parse_args()
	#session = requests.session()
	patentlist =  pd.read_csv(args.patent,header=0)['patent']
	## 模拟浏览器
	session.headers = { 
        "User-Agent": "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/106.0.0.0 Safari/537.36 "
    }
	## 企业版登录网址
	login_url = "https://db-server/external/authenticate/authorizations"
	login_data = { 
    'username': '15801517338',
    'password': 'mega123',
    'NVCVal' : "123"
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
	for patent in patentlist:
		Scrapy(session,patent,args.outdir)
		time.sleep(random.randint(2,5))
	loginout_url = "https://db-server/external/authenticate/authorization"
	session.delete(loginout_url)
	print ("Done")
	
def Scrapy(session,patent,outdir):
	#print (session.headers)
	print (patent)
	## 获取搜索的基本信息
	basic_url = "https://db-server/external/medical/databases/48/documents/" + patent
	basic_info = session.get(basic_url)
	basic_dict = json.loads(basic_info.text)
	response_result = pd.DataFrame(basic_dict)
	file = outdir + "/" + patent + ".txt"
	response_result.to_csv(file,sep="\t",encoding='utf8')
	#print (basic)

if __name__ == '__main__':
	main()
