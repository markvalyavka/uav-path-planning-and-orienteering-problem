import yaml

txt_filename = "/Users/markv/Desktop/Tsiligirides_2/tsiligirides_problem_2_budget_15.txt"
yaml_filename = "../input_configs/cfg_ts2_paper_benchmark.yaml"

start_location = end_location = None
locations = []
rewards = []

with open(txt_filename, "r") as f:
    # f.readline()
    # f.readline()
    tmax = float(f.readline().strip().split("\t")[0])
    lines = f.readlines()
    for l in lines:
        loc = list(map(float, l.strip().split("\t")))
        rewards.append(loc[2])
        loc[2] = 0
        locations.append(loc)
    start_location = locations[0]
    end_location = locations[1]
    locations.pop(0)
    locations.pop(0)
    rewards.pop(0)
    rewards.pop(0)

with open(yaml_filename, "r") as fa:
    yaml_cfg = yaml.full_load(fa)
    yaml_cfg['t_max'] = tmax
    yaml_cfg['start']['position'] = start_location
    yaml_cfg['end']['position'] = end_location
    yaml_cfg['locations'] = locations
    yaml_cfg['rewards'] = rewards

with open(yaml_filename, "w") as fa:
    yaml.dump(yaml_cfg, fa, default_flow_style=None)




