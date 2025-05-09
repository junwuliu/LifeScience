## 下面为借鉴和改造CFD得分程序
import os
import pickle
import pandas as pd
import numpy as np
import re
import sys
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)
def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open('{}/mismatch_score.pkl'.format(script_dir),'rb'))
        pam_scores = pickle.load(open('{}/pam_scores.pkl'.format(script_dir),'rb'))
        return (mm_scores,pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")
def calc_cfd(wt,sg,pam):
    mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)
