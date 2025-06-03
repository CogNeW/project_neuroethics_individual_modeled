# fatigue_HDDM_analysis.py - run hddm for task data from 'create_task_dataset.R
#
# Copyright (C) 2022 Kianté Fernandez, <kiantefernan@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 07/09/22      Kianté  Fernandez                   wrote initial code


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import hddm
import pickle

print(hddm.__version__)

data = hddm.load_csv("~/Dropbox/Medaglia_RiskBenefit_Neuroethics/DATA_ANALYSIS/data/HDDM_data.csv")

#include=('z')
m_temp0 = hddm.HDDM(data, p_outlier=.05, include=('z'))
m_temp1 = hddm.HDDM(data,  depends_on={'v': 'domain'}, p_outlier=.05,include=('z'))
m_temp2 = hddm.HDDM(data,  depends_on={'a': 'domain'}, p_outlier=.05, include=('z'))
m_temp3 = hddm.HDDM(data,  depends_on={'v': ['domain'], 'a': ['domain']}, p_outlier=.05, include=('z'))

m_temp0.find_starting_values()
m_temp1.find_starting_values()
m_temp2.find_starting_values()
m_temp3.find_starting_values()

m_temp0.sample(10000, burn=5000)
m_temp1.sample(10000, burn=5000)
m_temp2.sample(10000, burn=5000)
m_temp3.sample(10000, burn=5000)

m_temp0.dic
m_temp1.dic
m_temp2.dic
m_temp3.dic

pre_res0 = m_temp0.gen_stats()
pre_res0.to_csv('m0_results.csv')
pre_res1 = m_temp1.gen_stats()
pre_res1.to_csv('m1_results.csv')
pre_res2 = m_temp2.gen_stats()
pre_res2.to_csv('m2_results.csv')
pre_res3 = m_temp3.gen_stats()
pre_res3.to_csv('m3_results.csv')

#we need to code for a starting point bias in the model. 
#we know people are bias towards not selecting a treatment option
# m_risk_subj = hddm.HDDMRegressor(data, "v ~ risk + gain", include = 'z')
# 
# m_risk_subj.find_starting_values()
# m_risk_subj.sample(5000, burn=1000)
# m_risk_subj.dic
# pre_reg = m_risk_subj.gen_stats()
# 
# dont do this. Information is not supposed to effect boundaies
# m_risk_subj2 = hddm.HDDMRegressor(data, "a ~ risk + gain")
# m_risk_subj2.find_starting_values()
# m_risk_subj2.sample(5000, burn=1000)
# m_risk_subj2.dic
# pre_reg = m_risk_subj2.gen_stats()
#the bounds model fits worse than the drift regression
