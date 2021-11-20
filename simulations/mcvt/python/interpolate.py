import sys, subprocess
import glob
import numpy as np

from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton


# INTERPOLATE
def interpolate(list):
    print(list)
    x = [index for index, v in enumerate(np.ones(len(list)))]
    y = list
    f = interp1d(x, y)
    f2 = interp1d(x, y, kind='cubic')
    xnew = [index/100 for index, v in enumerate(np.ones(1+100*(len(list)-1)))]

    plt.plot(xnew, f2(xnew), '-')
    plt.legend(['Bottom', 'Iron', 'Top'], loc='best')
    plt.xlabel('Angle', fontsize=18)
    plt.ylabel('Torge Nm', fontsize=16)
    ax = plt.gca()
    ax.set_xlim([0, 45])
    ax.set_ylim([-20, 20])


## LOAD FILE
def loadFile(filename):
    f = open(filename, "r")
    data = f.read()
    f.close()

    lines = data.strip().split('\n')
    array = np.array([line.split(',') for line in lines])

    return np.array([[float(x) for x in lines] for lines in array.transpose()])

## NAVIGATION
counter = 0
def draw(index):
    global results, figure
    plt.cla()
    plt.clf()
    print(index)
    [interpolate(part) if index < 3 else None for index, part in enumerate(results[index])]
    figure.suptitle('Torgue/Angle for Iron angle ' + str(index), fontsize=20)

    figure.canvas.draw()

def onPrevious():
    global results, counter
    print("P")
    counter = counter - 1 if counter > 0 else counter
    draw(counter)

def onNext():
    global results, counter
    print("N")
    counter = (counter + 1) if counter < len(results) - 1 else counter
    draw(counter)
def onButton(event):
    if event.button == MouseButton.LEFT:
        onNext()
    elif event.button == MouseButton.RIGHT:
        onPrevious()

figure = plt.figure()
figure.canvas.mpl_connect('button_press_event', onButton)

## MAIN CODE

if len(sys.argv) < 2:
    print("Missing parameter filename!")
    exit(1)

filepaths = sys.argv[1:]
print("Loading " + str(len(filepaths)) + " files")
results = [loadFile(path) for path in filepaths]

draw(counter)
plt.show()


