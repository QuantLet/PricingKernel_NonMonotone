# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:07:59 2024

@author: Ruting Wang


"""
from PIL import Image

import os
import pandas as pd
# import datetime
# import statsmodels.api as sm
# import matplotlib.pyplot as plt
import numpy as np 
# import time
# import sys
# import pdb
# import pickle
# import gc
# import csv
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil
from datetime import datetime

path = "/Users/ruting/Library/Mobile Documents/com~apple~CloudDocs/PK_BTC/Upload_quantlet/PricingKernel_NonMonotone/CDIPlot"
os.chdir(path)

# 打开图片
png_file = "BRC_EPK_13day_2022-10-29"
image = Image.open(png_file+".png").convert("RGBA")

# 获取图片的像素数据
data = image.getdata()

# 替换白色背景为透明
new_data = []
for item in data:
    # 如果是白色像素 (255, 255, 255)，设置为透明
    if item[:3] == (255, 255, 255):  # RGB 值
        new_data.append((255, 255, 255, 0))  # 设置 alpha 通道为 0
    else:
        new_data.append(item)

# 更新图片数据
image.putdata(new_data)

# 保存为带透明背景的新图片
image.save(png_file+"_transparent.png", "PNG")
