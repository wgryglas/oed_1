__author__ = 'wgryglas'
# -*- coding: utf-8 -*-


def perform(dirs, files, par):
    """
    Script which imports data from excel sheet into
    numpy data files for faster access.
    """
    import pyexcel as pe
    import numpy as np
    from wg.tools.system import ask_to_clean_dir
    import os, sys, re

    reload(sys)  # Reload does the trick!
    sys.setdefaultencoding('UTF8')


    ask_to_clean_dir(dirs.expriment)


    sheet = pe.get_sheet(file_name=files.experiment_spreadsheet, sheet_name=par.experiment_sheet_name)

    infos = {sheet['F1']: sheet['G1'], sheet['F2']: sheet['G2'], sheet['F3']: sheet['G3'],
             sheet['H1']: sheet['I1'], sheet['H2']: sheet['I2'], sheet['H3']: sheet['I3'],
             sheet['J1']: sheet['K1'], sheet['F4']: sheet['G4']}

    np.save(dirs.experiment+os.sep+'infos', infos)

    channels = sheet.column[4][7:]
    np.save(dirs.experiment+os.sep+'channels', channels)

    pnts = np.array([sheet.column[2][7:], sheet.column[3][7:]]).T
    np.save(dirs.experiment+os.sep+'points', pnts)


    rPressure = re.compile(r'([\d\.]+)')

    for c in range(5, sheet.number_of_columns(), 4):
        pressure = re.search(rPressure, sheet.column[c][4]).group(0).replace('.', '-')
        outname = 'press'+pressure+par.experiment_output_suffix
        data = np.array([sheet.column[c+i][7:] for i in range(4)]).T
        np.save(dirs.experiment+os.sep+outname, data)