import matplotlib.pyplot as plt

def graph(x,y):
    
    #Initialize the plot
    fig = plt.figure()

    #Set up axes
    ax = fig.add_subplot(111)

    #Activate grid for the plot
    ax.grid(True)

    #Set the names to the axes
    ax.set(title='k vs |f(x)|', xlabel='k', ylabel='|f(x)|')

    #Plot the data
    ax.plot(x, y, color='blue', linewidth=1.5)

    #Show the plot
    plt.show()
