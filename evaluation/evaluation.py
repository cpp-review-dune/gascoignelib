import os
import sys
import unittest
import subprocess as sp
import os
import json

class dir:
    def __init__(self, root):
        if root[-1] == "/":
            root = root[:-1]
        self.dir = root.split("/")
    
    def push(self, folder):
        self.dir.append(folder)
    
    def pop(self):
        self.dir.pop()
    
    def __str__(self):
        return '/'.join(self.dir) + "/"

def insert(list, key, value):
    if len(key) > 1:
        if not key[0] in list:
            list[key[0]] = {}
        insert(list[key[0]], key[1:], value)
    else:
        list[key[0]] = value

def inData(list, key):
    if len(key) == 0:
        return True
    if key[0] in list:
        return inData(list[key[0]], key[1:])
    return False

def get(list, key):
    if len(key) > 1:
        return get(list[key[0]], key[1:])
    return list[key[0]]
    
def runExamples():
    root_dir = dir(os.getcwd())
    evaluation_dir = dir(str(root_dir) + "evaluation/")

    commit = sp.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True)
    builddir = dir(str(root_dir) + "build/")
    outputfolder = dir(str(builddir) + "eval-" + (commit.stdout).strip() + "/")


    if not os.path.isdir(str(builddir)):
        os.mkdir(str(builddir))
    else:
        # if os.path.isdir(str(outputfolder)):
        #     return str(outputfolder)
        sp.call("rm -r build/*", shell=True)
        
    os.mkdir(str(outputfolder))

    examples = {}
    with open(str(evaluation_dir) + 'examples.json', 'r') as f:
        examples = json.load(f)
    sp.run(["cp", str(evaluation_dir) + 'examples.json', str(outputfolder) + "exmaples.json"])

    os.chdir('build')

    for device in ["CPU", "GPU"]:
    # for device in [ "GPU", ]:
        print("Building " + device)
        os.mkdir(device)
        os.chdir(device)
        with open(str(outputfolder) + 'build_' + device + '.out', 'w') as f:
            sp.run(["date"], stdout=f)
            sp.run(["git", "log", "-n", "1"], stdout=f)
            if "GPU" in device:
                sp.run(["cmake", '-DUSE_CUDA=ON', '../..'], stdout=f)
            else:
                sp.run(["export", 'OMP_NUM_THREADS=8'], shell=True, stdout=f)
                sp.run(["cmake", '-DUSE_CUDA=OFF', '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -fopenmp', '-DCMAKE_EXE_LINKER_FLAGS_RELEASE=-fopenmp', '../..'], stdout=f)
            sp.run(["make", '-j'], stdout=f)

        os.chdir('evaluation/examples')
        for folder, example in examples.items():
            print("Running: " + folder)
            outputfolder.push(folder)
            if not os.path.isdir(str(outputfolder)):
                os.mkdir(str(outputfolder))
            os.chdir(folder)
            for file in example["copyfiles"]:
                sp.run(["cp", str(root_dir) + 'evaluation/examples/' + folder + "/" + file, "."])
            with open(str(evaluation_dir) + example["parameterfile"], 'r') as f:
                paramfile = f.read()
            for preref in range(1,example["prerefine"]):
                print("Prerefine-level: " + str(preref))
                paramfile_out = paramfile.replace("{{prerefine}}", str(preref))
                with open(example["parameterfile"], 'w') as f:
                    f.write(paramfile_out)
                with open(str(outputfolder) + "/" + folder + "_" + device + "_" + str(preref) + ".txt", "w") as f:
                    f.write("*** Parafile ***\n")
                    f.write(paramfile_out)
                    f.write("\n*** Output ***\n")
                    f.flush()
                    sp.run("../../../bin/" + example["executable"], stdout=f, stderr=f)
            os.chdir("../")
            outputfolder.pop()
        os.chdir("../../../")
    os.chdir("../")
    return str(outputfolder)

def read(root):
    data = {}
    for file in os.listdir(root):
        path = root + "/" + file
        if os.path.isfile(path):
            pass
        elif os.path.isdir(path):
            for parafile in os.listdir(path):
                experiment = parafile.split(".")[0].split("_")
                ignore = True
                with open(path + "/" + parafile, "r") as f:
                    for line in f:
                        line = line.strip()
                        if "*** Timer ***" in line:
                            ignore = False
                            continue
                        if ignore:
                            continue
                        timer = line[0:26].strip()
                        times = [a for a in line[26:].split(" ") if a]
                        if len(times) > 0:
                            if len(times) == 4:
                                insert(data, [experiment[0], timer.strip(), int(experiment[2]), experiment[1], "time"], float(times[0]))
                            insert(data, [experiment[0], timer.strip(), int(experiment[2]), experiment[1], "count"], times[-1])
    return data

def generateKeys(data):
    keys = []
    for examples, examples_data in data.items():
        for metrics, metric_data in examples_data.items():
            for refinement, ref_data in metric_data.items():
                for device, device_data in ref_data.items():
                    for type, type_data in device_data.items():
                        keys.append([examples,metrics,refinement, device, type])
    return keys

class TestGascogine(unittest.TestCase):

    def setUp(self):
        self.eps = 2
        self.reference = read(os.getcwd() + "/evaluation/reference/")
        self.keys = generateKeys(self.reference)
        self.test_folder = runExamples()
        self.test_data = read(self.test_folder)

    def test_metrics(self):
        for key in self.keys:
            with self.subTest(key=key):
                self.assertTrue(inData(self.test_data, key))
                if "count" in key[-1]:
                    self.assertEqual(get(self.reference, key), get(self.test_data, key))
                if "time" in key[-1]:
                    self.assertTrue( (get(self.test_data, key) / get(self.reference, key)) < 1 + self.eps)


if __name__ == '__main__':
    root = sys.argv[1]
    data = TestGascogine.read(root)
    print(str(data))
    runner = unittest.main()

