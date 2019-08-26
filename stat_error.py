import os
import sys
import re

class TCA:
    def __init__(self, benchmark, tuple_size, tabu_step_size=None, random_size=None):
        self.benchmark = benchmark
        self.tuple_size = tuple_size
        self.tabu_step_size = tabu_step_size
        self.random_size = random_size
        self.run_times = []
        self.last_sizes = []

    #def __hash__(self):
    #    if self.tabu_step_size is not None:
    #        return hash("_".join([self.benchmark, str(self.tuple_size), str(self.tabu_step_size), str(self.random_size)]))
    #    else:
    #        return hash("_".join([self.benchmark, str(self.tuple_size)]))

    def get_key(self):
        if self.tabu_step_size is not None:
            return "_".join([self.benchmark, str(self.tuple_size), str(self.tabu_step_size), str(self.random_size)])
        else:
            return "_".join([self.benchmark, str(self.tuple_size)])

    def get_seq(self):
        return self.benchmark.split("_")[1]

    def __cmp__(self, s):
        this_seq = self.get_seq()
        s_seq = s.get_seq()
        if s_seq.isdigit() and this_seq.isdigit():
            return cmp(int(this_seq), int(s_seq))
        elif this_seq.isdigit():
            return 1
        elif s_seq.isdigit():
            return -1
        else:
            return cmp(this_seq, s_seq) 


    #def __eq__(self, other):
    #    if not ( self.benchmark == other.benchmark and self.tuple_size == other.tuple_size):
    #        return False
    #    if self.tabu_step_size is not None and not (self.tabu_step_size == other.tabu_step_size and self.random_size == other.random_size):
    #        return False
    #    return True

    @staticmethod
    def build_tca(filename):
        str_arr = filename.split("_")
        if len(str_arr) != 4 and len(str_arr) != 6:
            return None

        name = str_arr[0] + "_" + str_arr[1]
        name = name.replace("error.", "")
        seed = int(str_arr[2])
        tuple_size = 2 if str_arr[-1] == "2way" else 3
        if len(str_arr) == 6:
            tabu_step_size = int(str_arr[3])
            random_size = int(str_arr[4])
            return TCA(name, tuple_size, tabu_step_size, random_size)
        else:
            return TCA(name, tuple_size)

    def add_tca_value(self, t, s):
        self.last_sizes.append(s)
        self.run_times.append(t)

    def print_info(self):
        total_size = 0 
        best_size = sys.maxint
        total_time = 0.0
        for s in self.last_sizes:
            total_size += s
            if best_size > s:
                best_size = s
        for t in self.run_times:
            total_time += t
        if len(self.last_sizes) ==0:
            avg_size = avg_time = 0
        else:
            avg_size = total_size*1.0/len(self.last_sizes)
            avg_time = total_time*1.0/len(self.last_sizes)
        if self.random_size is None:
            return ",".join([self.benchmark, str(avg_time), str(avg_size), str(best_size)])
        else:
            return ",".join([self.benchmark, str(self.tabu_step_size), str(self.random_size),str(avg_time), str(avg_size), str(best_size)])

def stat_file(tca_dict, filename):
    tca = TCA.build_tca(filename)
    if tca is None:
        return False
    key = tca.get_key() 
    if key not in tca_dict:
        tca_dict[key] = tca
    t = None
    s = None
    with open(filename, "r") as fh:
        last_line = ""
        for line in fh.readlines():
            line = line.strip()
            if line.find("CoveringArray::optimize") > 0:
                t, s = re.split("\s+", last_line)
                #t, s = last_line.split(" ")
                #print name, seed, tuple_size, t , s
                break
            if line.find("end") < 0:
                last_line = line
    print filename, t, s
    tca_dict[key].add_tca_value(float(t),int(s))
    #return [name, seed, tuple_size, t, s] 
    return True

tca_dict = {}
for filename in os.listdir(os.getcwd()):
    if not filename.startswith("error"):
        continue
    stat_file(tca_dict, filename)

if len(tca_dict) == 0:
    sys.exit()
l = tca_dict.values()
l.sort()
fh_2way = open("result_2way.csv", "w")
fh_3way = open("result_3way.csv", "w")
if l[0].random_size is not None:
    fh_2way.write("benchmark,tabu_step_size,random_size,avg_time,avg_size,best_size\n")
    fh_3way.write("benchmark,tabu_step_size,random_size,avg_time,avg_size,best_size\n")
else:
    fh_2way.write("benchmark,avg_time,avg_size,best_size\n")
    fh_3way.write("benchmark,avg_time,avg_size,best_size\n")
for tca in l:
    info = tca.print_info()
    if tca.tuple_size == 2 :
        fh_2way.write(info+"\n")
    elif tca.tuple_size == 3:
        fh_3way.write(info+"\n")
fh_2way.close()
fh_3way.close()
