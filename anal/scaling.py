import timePerf
import argparse
import matplotlib.pylab as plt
import numpy as np

def concatContentFiles(fileList):

    output=""
    for file in fileList:
        with open(file) as f:
            output+=f.read()
    return output

def plotScaling(datas):

  

    for label,data in datas.groupby("label"):
        x=data["cores"]

        ncMin=np.min(data["cores"])
        print (ncMin)
        t0=np.max(data[ data["cores"] == ncMin ]["time"]) 
        #t0=np.max(data["time"])
        y=t0/data["time"]


        plt.plot(x,y,"o",label=label)

        plt.plot(x,x/ncMin,"--")

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create scaling plot')
    parser.add_argument('filenames',nargs="+", help='List of files to parse') 

    args  = parser.parse_args()
    data= timePerf.generateDataFromStringOutput(  concatContentFiles( args.filenames) )
    print(data)
    plotScaling(data)

    plt.legend()
    plt.show() 