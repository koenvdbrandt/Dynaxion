import sys
import os
import random
import numpy as np

from array import array



eventArr = array('i', [0])
energyArr = array('d', [0])
timeArr = array('d', [0])
positionxArr = array('d', [0])
positionyArr = array('d', [0])
positionzArr = array('d', [0])
pdg_codeArr = array('i', [0])
track_idArr = array('i', [0])
parent_idArr = array('i', [0])


# Class to store deposited charges in
class depositedCharge:
    
    def __init__(self):
        pass

    def setEventNr(self, eventNr):
        self.eventNr = eventNr
    def setEnergy(self, energy):
        self.energy = energy
    def setTime(self, time):
        self.time = time
    def setPositionx(self, positionx):
        self.positionx = positionx
    def setPositiony(self, positiony):
        self.positiony = positiony
    def setPositionz(self, positionz):
        self.positionz = positionz
    def setDetector(self, detector):
        self.detector = detector
    def setPdg_Code(self, pdg_code):
        self.pdg_code = pdg_code
    def setTrack_Id(self, track_id):
        self.track_id = track_id
    def setParent_Id(self, parent_id):
        self.parent_id = parent_id
    

    # Required for TTrees
    def fillDepositionArrays(self):

        global eventArr, energyArr, timeArr, positionxArr, positionyArr, positionzArr, pdg_codeArr, track_idArr, parent_idArr
    
        eventArr[0] = self.eventNr
        energyArr[0] = self.energy
        timeArr[0] = self.time
        positionxArr[0] = self.positionx
        positionyArr[0] = self.positiony
        positionzArr[0] = self.positionz
        pdg_codeArr[0] = self.pdg_code
        track_idArr[0] = self.track_id
        parent_idArr[0] = self.parent_id


    # Required for CSV output
    def getDepositionText(self):

        text = str(self.pdg_code) + ", "
        text += str(self.time) + ", "
        text += str(self.energy) + ", "
        text += str(self.positionx) + ", "
        text += str(self.positiony) + ", "
        text += str(self.positionz) + ", "
        text += str(self.detector) + ", "
        text += str(self.track_id) + ", "
        text += str(self.parent_id)

        text += "\n"

        return text
        

# Calculation of straight particle trajectories in the sensor
def createParticle(nsteps):

    pdgCode = 11

    # For time calculation in ns
    timeOffset = 1.

    # Arbitrary entry and exit points in global coordinates
    entryX = random.gauss(2., 0.5)
    exitX = random.gauss(entryX, 0.05)
    entryY = random.gauss(-0.5, 0.5)
    exitY = random.gauss(entryY, 0.05)
    entryPoint = np.array([entryX, entryY, -0.1420])
    exitPoint = np.array([exitX, exitY, 0.1420])

    distVect = exitPoint - entryPoint
    totalDist = np.sqrt(distVect.dot(distVect))

    # Calculate mean energy deposition
    eVPermm = 390/0.001
    eVPerStep = eVPermm * (totalDist/nsteps)

    # Track parametrization
    zPositions = np.linspace(entryPoint[2],exitPoint[2],nsteps, endpoint=True)
    slope = np.array([(exitPoint[0]-entryPoint[0])/(exitPoint[2]-entryPoint[2]), (exitPoint[1]-entryPoint[1])/(exitPoint[2]-entryPoint[2])])
    offset = np.array([entryPoint[0]-slope[0]*entryPoint[2], entryPoint[1]-slope[1]*entryPoint[2]])

    # Create vector of deposited charges
    deposits = []
    
    for zPos in zPositions:
        xyPos = offset + slope*zPos
        xyzPos = np.append(xyPos,zPos)

        # Calculate time in ns
        time = timeOffset + zPos*1e-3 / 3e8 * 1e9
        energy = random.gauss(eVPerStep,eVPerStep/30.)

        # Create deposited charge as an object
        deposit = depositedCharge()
        deposit.setEnergy(energy)
        deposit.setTime(time)
        deposit.setPositionx(xyzPos[0])
        deposit.setPositiony(xyzPos[1])
        deposit.setPositionz(xyzPos[2])
        deposit.setPdg_Code(pdgCode)
        # Track and parent ID arbitrary in this case
        deposit.setTrack_Id(0)
        deposit.setParent_Id(0)

        # Push back vector of deposited charges
        deposits.append(deposit)

    return deposits
    
def user_input(question):
    if sys.version_info.major == 3:
        return input(question)
    elif sys.version_info.major == 2:
        return raw_input(question)
    else:
        print("Python version could not be determined.")
        exit(1)
    


if __name__ == '__main__':

    # Check for availability of ROOT
    rootAvailable = True
    try:
        import ROOT
    except:
        print("ROOT unavailable. Install ROOT with python option if you are interested in writing TTrees.")
        rootAvailable = False

    
    # Ask whether to use TTrees or CSV files
    if rootAvailable:
        writeOption = user_input("Generate TTrees (a), a CSV file (b) or both (c)? ")
        writeROOT = False
        writeCSV = False
        if writeOption=="a":
            writeROOT = True
        elif writeOption=="b":
            writeCSV = True
        elif writeOption=="c":
            writeROOT = True
            writeCSV = True
        else:
            print("Use one of the three options.")
            exit(1)
    else:
        print("Will just write a CSV file.")
        writeCSV = True
        writeROOT = False
    

    filenamePrefix = "deposition"
    
    rootfilename = filenamePrefix + ".root"
    csvFilename = filenamePrefix + ".csv"
    
    # Define detector name
    detectorName = user_input("Name of your detector: ")

    # Ask for the number of events
    events = int(user_input("Number of events to process: "))

    # Ask for the number of steps along the track:
    nsteps = int(user_input("Number of steps along the track in the sensor: "))

    if writeROOT:
        # Open the file and create the tree
        rootfile = ROOT.TFile(rootfilename,"RECREATE")
        tree = ROOT.TTree("treeName","treeTitle")

        detectorArr = ROOT.vector('string')()
        detectorArr.push_back(detectorName)

        # Create the branches
        eventBranch = tree.Branch("event", eventArr, "event/I")
        energyBranch = tree.Branch("energy", energyArr, "energy/D")
        timeBranch = tree.Branch("time", timeArr, "time/D")
        positionxBranch = tree.Branch("position.x", positionxArr, "position.x/D")
        positionyBranch = tree.Branch("position.y", positionyArr, "position.y/D")
        positionzBranch = tree.Branch("position.z", positionzArr, "position.z/D")
        # The char array for the detector branch is a bit more tricky, since you have to give it the length of the name
        detectorBranch = tree.Branch("detector", detectorArr[0], "detector["+str(len(detectorArr[0]))+"]/C")
        pdg_codeBranch = tree.Branch("pdg_code", pdg_codeArr, "pdg_code/I")
        track_idBranch = tree.Branch("track_id", track_idArr, "track_id/I")
        parent_idBranch = tree.Branch("parent_id", parent_idArr, "parent_id/I")
        

    if writeCSV:
        fout = open(csvFilename,'w')
    

    for eventNr in range(0,events):
        print("Processing event " + str(eventNr))

        # Get the depositions for the particle created
        deposits = createParticle(nsteps)

        if writeCSV:
            # Write the event number
            text = "\nEvent: " + str(eventNr) + "\n"
            fout.write(text)

        for deposit in deposits:
            # Add information to the depositions
            deposit.setEventNr(eventNr)
            deposit.setDetector(detectorName)

            if writeROOT:
                # Fill the arrays and then write them to the tree
                deposit.fillDepositionArrays()
                tree.Fill()

            if writeCSV:
                # extract the text lines for the individual depositions
                text = deposit.getDepositionText()
                fout.write(text)


    if writeROOT:
        # Inspect tree and write ROOT file
        tree.Scan()
        rootfile.Write()
        rootfile.Close()

    if writeCSV:
        # End the file with a line break to prevent from the last line being ignored due to the EOF
        fout.write("\n")
        fout.close()

        
