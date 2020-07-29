# -*- coding: utf-8 -*-
# @Author: Wangj
# @Date:   2019-12-05 17:18:59
# @Last Modified by:   Hanfeng
# @Last Modified time: 2019-12-05 17:32:08

import os
import time
basedir = os.getcwd()
filelists = []
# 指定想要统计的文件类型
whitelist = ['cpp', 'hpp']
#遍历文件, 递归遍历文件夹中的所有


def getFile(basedir):
    global filelists
    for parent, dirnames, filenames in os.walk(basedir):
        for filename in filenames:
            ext = filename.split('.')[-1]
            #只统计指定的文件类型，略过一些log和cache文件
            if ext in whitelist:
                fullDir = os.path.join(parent, filename)
                relPath = os.path.relpath(fullDir, basedir)
                filelists.append(relPath)


#统计一个文件的行数
def countLine(fname):
    count = 0
    for file_line in open(fname).readlines():
        if file_line != '' and file_line != '\n':  #过滤掉空行
            count += 1
    print(fname + '----', count)
    return count


if __name__ == '__main__':
    startTime = time.clock()
    getFile(basedir)
    totalline = 0
    for filelist in filelists:
        totalline = totalline + countLine(filelist)
    print('total lines:', totalline)
    print('Done! Cost Time: %0.2f second' % (time.clock() - startTime))
