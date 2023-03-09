import subprocess
from subprocess import Popen
from subprocess import Popen,PIPE,STDOUT


t_maxes = ["10", "15", "20", "25", "30", "35", "40"]
t_maxes = ["10", "15", "20"]
Vs = ["1", "3", "6"]

results = []

for t_max in t_maxes:
     for V in Vs:

          # heu = Popen(["../main", t_max, V], stdout=subprocess.PIPE)
          results.append(Popen(["../main", t_max, V], stdout=subprocess.PIPE))
          # res = heu.stdout.read().decode('utf-8')
          # print(res[res.find("<result>")+8:res.find("</result>")])

for r in results:
     res = r.stdout.read().decode('utf-8')
     print(res[res.find("<result>") + 8:res.find("</result>")])