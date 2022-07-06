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
    
def runTest():
    root_dir = dir(os.getcwd())
    evaluation_dir = dir(str(root_dir) + "evaluation/")

    commit = sp.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True)
    outputfolder = dir(str(root_dir) + "build/eval-" + (commit.stdout).strip() + "/")

    examples = {}
    with open(str(evaluation_dir) + 'examples.json', 'r') as f:
        examples = json.load(f)

    sp.call("rm -r build", shell=True)
    os.mkdir('build')
    os.chdir('build')

    if not os.path.isdir(str(outputfolder)):
        os.mkdir(str(outputfolder))

    for device in ["CPU", "GPU"]:
        print("Building " + device)
        with open(str(outputfolder) + 'build_' + device + '.out', 'w') as f:
            sp.run(["git", "log", "-n", "1"], stdout=f)
            if "GPU" in device:
                sp.run(["cmake", '-DUSE_CUDA=ON', '..'], stdout=f)
            else:
                sp.run(["export", 'OMP_NUM_THREADS=8'], shell=True, stdout=f)
                sp.run(["cmake", '-DUSE_CUDA=OFF', '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -fopenmp', '-DCMAKE_EXE_LINKER_FLAGS_RELEASE=-fopenmp', '..'], stdout=f)
            sp.run(["make", '-j'], stdout=f)


        os.chdir('examples')
        for dir, example in examples.items():
            print("Running: " + dir)
            outputfolder.push(dir)
            if not os.path.isdir(str(outputfolder)):
                os.mkdir(str(outputfolder))
            os.chdir(dir)
            for file in example["copyfiles"]:
                sp.run(["cp", str(root_dir) + 'examples/' + dir + "/" + file, "."])
            with open(str(evaluation_dir) + example["parameterfile"], 'r') as f:
                paramfile = f.read()
            for preref in range(1,example["prerefine"]):
                print("Prerefine-level: " + str(preref))
                paramfile_out = paramfile.replace("{{prerefine}}", str(preref))
                with open(example["parameterfile"], 'w') as f:
                    f.write(paramfile_out)
                with open(str(outputfolder) + "/" + dir + "_" + device + "_" + str(preref) + ".txt", "w") as f:
                    f.write("*** Parafile ***\n")
                    f.write(paramfile_out)
                    f.write("\n*** Output ***\n")
                    f.flush()
                    sp.run("./" + example["executable"], stdout=f, stderr=f)
            os.chdir("../")
            outputfolder.pop()
        os.chdir("../")
    return str(outputfolder)
