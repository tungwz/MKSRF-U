#!/usr/bin/python

import numpy
import os

#year = 2000
reg = numpy.loadtxt('reg_5x5')
cnt = 1


for i in reg:
    print(f'processing {cnt}, {i[0]}, {i[1]}, {i[2]}, {i[3]}')

    #print('ncks -h'
    #        + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
    #        + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
    #        + ' CoLM_UCPs_500m.nc -o /tera10/yuanhua/dongwz/urban_data/CoLMUCPs_5x5_v2/RG_'
    #        + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
    #        + '.UCPs' + '.nc')

    os.system('ncks -h'
            + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
            + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
            + ' CoLM_UCPs_ROOF.nc -O /tera12/yuanhua/data/CoLMrawdata/urban_morphology/roof_height_fraction_Liao/RG_'
            + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
            + '.ROOF1km.Liao' + '.nc')

    os.system('ncks -h'
            + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
            + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
            + ' CoLM_UCPs_HL.nc -O /tera12/yuanhua/data/CoLMrawdata/urban_morphology/building_HL_Liao/RG_'
            + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
            + '.HL1km_Liao' + '.nc')

    os.system('ncks -h'
            + ' -d lat,' + str(int(i[2])) + '.,' + str(int(i[0])) + '.'
            + ' -d lon,' + str(int(i[1])) + '.,' + str(int(i[3])) + '.'
            + ' CoLM_UCPs_Fgper.nc -O /tera12/yuanhua/data/CoLMrawdata/urban_morphology/pervious_fraction_Liao/RG_'
            + str(int(i[0])) + '_' + str(int(i[1])) + '_' + str(int(i[2])) + '_' + str(int(i[3]))
            + '.Fgper1km_Liao' + '.nc')

    cnt = cnt + 1
