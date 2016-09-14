__author__ = 'wgryglas'

from wg.tools.function import assert_arg_iterable
from exceptions import IndexError
import re
import numpy as np


class _DataStore:
    def __init__(self):
        self.name = None
        self.data = None
        self.size = -1
        self.id = 0


class RedBoundaryReader:
    def __init__(self, fname):
        print "reading"+fname
        self.fname = fname
        self.data = {}
        self.__vars = None
        self.__read_fun = self.__read_title
        self.__curr_data = _DataStore()
        self.surfaces = list()

    def read(self):
        name=""
        with open(self.fname, "r") as f:
            for l in f:
                self.__read_fun(l)

        return self

    _reTitle = re.compile("TITLE\s*=\s*\"(.+)\"")

    def __read_title(self, line):
        res = self._reTitle.search(line)
        if res:
            self.__curr_data.name = res.group(1)
            self.__read_fun = self.__read_variables

    _re_variables_find = re.compile('VARIABLES\s*=\s*(.*)')
    _re_variables_split = re.compile("\"([^\"]*)\"")

    def __read_variables(self, line):
        if self.__vars is None:
            self.__vars = self._re_variables_split.findall(self._re_variables_find.search(line).group(1))
        self.__read_fun = self.__read_size

    _re_size = re.compile("N\s*=\s*([0-9]+)")

    def __read_size(self, line):
        self.__curr_data.size = int(self._re_size.findall(line)[0])
        self.__curr_data.data = np.zeros((self.__curr_data.size, len(self.__vars)))
        self.__read_fun = self.__read_data

    def __read_data(self, line):
        if self.__curr_data.id >= self.__curr_data.size:
            self.__curr_data.id = 0
            self.__read_fun = self.__read_edges
            self.data[self.__curr_data.name] = self.__curr_data.data
            self.surfaces.append(self.__curr_data.name)
            self.__curr_data = _DataStore()
        else:
            self.__curr_data.data[self.__curr_data.id] = np.fromstring(line, sep=' ')
            self.__curr_data.id += 1

    def __read_edges(self, line):
        if self.__curr_data.id >= self.__curr_data.size:
            self.__read_fun = self.__read_title
        else:
            self.__curr_data.id += 1


    def write(self, fname, vars, surfaces=None):
        if not surfaces:
            surfaces = [self.surfaces[0]]

        varIds = self.variableIds(vars)

        with open(fname,'w') as f:
            title = ''
            for s in surfaces:
                if len(surfaces) > 1:
                    title += 'TITLE="'+s+'"\t'
                title += 'N=' + str(self.data[s].shape[0])+'\n'
                f.write(title)
                for row in self.data[s]:
                    f.write('\t'.join(map(str,row[varIds]))+'\n')

    @property
    def variables(self):
        return self.__vars

    @assert_arg_iterable(args=1, kwargs='names')
    def variableIds(self, names):
        ids = list()
        for i, v in enumerate(self.__vars):
            if v in names:
                ids.append(i)

        if len(ids) == 1:
            return ids[0]
        else:
            return ids

    def titles(self):
        return self.data.keys()

    def __getitem__(self, item):
        return self.data[item]

    def __iter__(self):
        return self.data.__iter__()


def write_coordinates(out_file, X, Y, sep='\t'):

    with open(out_file,'w') as f:
        f.write("N=%d" % len(X))
        for x, y in zip(X, Y):
            f.write(("\n%lf"+sep+"%lf") % (x, y))


def readBoundaryNodes(nodesFile, surfaceName=None, sep='\t'):

    data = {'N':0, 'nid':0, 'x':None, 'y':None, 'func':None}

    def readXY(l):
        xy = np.fromstring(l, sep=sep)
        data['x'][data['nid']] = xy[0]
        data['y'][data['nid']] = xy[1]
        data['nid'] += 1

        return data['nid'] < data['N']


    def readSize(l):
        N = int(re.findall('N=([0-9]+)',l)[0])
        data['N'] = N
        data['x'] = np.zeros(N)
        data['y'] = np.zeros(N)
        data['func'] = readXY
        return True

    data['func'] = readSize

    with open(nodesFile, 'r') as f:
        for l in f:
            if not data['func'](l):
                break

    return data['x'], data['y']

    # if not surfaceName:





