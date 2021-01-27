from math import *

INF = float("inf")
def log_sum_exp(numlist):
    # using log-trick to compute log(sum(exp(x) for x in numlist))
    # mitigate the problem of underflow
    try:
        minx = min([x for x in numlist if x != -INF])
    except:
        return -INF    
    s = sum(exp(x-minx) for x in numlist)
    return minx + log(s) if s > 0 else -INF

class HMM:
    def __init__(self):
        #the alphabet of the HMM model. This is a list of characters.
        self.alphabet = []
        # emission probability of the I states, one dictionary for each I in the model
        self.eI = [] 
        # emission probability of the M states, one dictionary for each M in the model
        # the first M state is called 'B'; it never emits anything so the associated dictionary is always empty
        self.eM = [{}] 
        # transition probability, one dictionary for each set of states (D,M,I) in the model
        self.t = [] 
        # number of matching states in the model, excluding B and E
        self.nstate = 0
    
    def load(self,hmmfile):
    # only load the first model in the given hmmfile if there are two or more models
        with open(hmmfile,'r') as fin:
            for line in fin:
                stream = line.strip().split()
                if stream[0] == "LENG":
                    self.nstate = int(stream[1])
                if stream[0] == "HMM": 
                    # read alphabet
                    self.alphabet = stream[1:]
                    # read transition order
                    stream = fin.readline().strip().split()
                    trans_order = [(y[0]+y[3]).upper() for y in stream]
                    # read the next line, if it is the COMPO line then ignore it and read one more line
                    stream = fin.readline().strip().split()
                    if stream[0] == "COMPO":
                        stream = fin.readline().strip().split()
                    # now the stream should be at the I0 state; read the emission of the I0 state 
                    e = {}
                    for (x,y) in zip(self.alphabet,stream):
                        e[x] = -float(y)
                    self.eI.append(e)    
                    # now the stream should be at the B state; read in the transition probability
                    stream = fin.readline().strip().split()
                    tB = {'MM':-INF,'MD':-INF,'MI':-INF,'IM':-INF,'II':-INF,'ID':-INF,'DM':-INF,'DI':-INF,'DD':-INF}
                    for x,y in zip(trans_order,stream):
                        tB[x] = -INF if y == '*' else -float(y)
                    self.t.append(tB)
                    break    
            
            for i in range(1,self.nstate+1):
                # read each set of three lines at a time
                stream = fin.readline().strip().split() # this one is the emission of the M state
                if float(stream[0]) != i:
                    print("Warning: incosistent state indexing in hmm file; expecting state "+ str(i) + "; getting state " + stream[0])
                e = {}
                for x,y in zip(self.alphabet,stream[1:]):
                    e[x] = -INF if y == "*" else -float(y) 
                self.eM.append(e)    
                # the next line is the emission of the I state
                stream = fin.readline().strip().split() 
                e = {}
                for x,y in zip(self.alphabet,stream):
                    e[x] = -INF if y == "*" else -float(y) 
                self.eI.append(e)                    
                # the next line contains the transition probs
                stream = fin.readline().strip().split() # this one is the transition prob
                tB = {'MM':-INF,'MD':-INF,'MI':-INF,'IM':-INF,'II':-INF,'ID':-INF,'DM':-INF,'DI':-INF,'DD':-INF}
                for x,y in zip(trans_order,stream):
                    tB[x] = -INF if y == '*' else -float(y)
                self.t.append(tB)
    
    def compute_llh(self,query):
        # compute likelihood of an aligned query to the HMM
        # return -inf if the query is not properly aligned
        j = 0
        prev_state = 'M'
        llh = 0

        for c in query:
            if c == '-': # gap -> deletion
                curr_state = 'D'
            elif c >= 'a' and c <= 'z': # lowercase -> insertion
                curr_state = 'I'
            elif c >= 'A' and c <= 'Z': # uppercase -> match
                curr_state = 'M'
            else: # encounter invalid symbol
                return -INF
            
            trans = prev_state + curr_state

            # penalize the transition
            if trans in self.t[j]:
                llh += self.t[j][trans]
            else: # encounter invalid transition
                return -INF 
            
            # transit: update the state index
            j += (curr_state != 'I') # move to the next state unless this is a move towards an 'I'
            if j > self.nstate: # reach end of the HMM chain but not end of the sequence
                return -INF

            # penalize the emission
            if curr_state == 'M':
                llh += self.eM[j][c]
            elif curr_state == 'I': 
                llh += self.eI[j][c.upper()]
            
            # update state
            prev_state = curr_state

        if j != self.nstate: # does not reach the E state at the end
            return float("-inf")
        
        # the last move: towards the 'E' state
        trans = prev_state + 'M'
        llh += self.t[j][prev_state+'M']    
        return llh      
    '''
    #the alphabet of the HMM model. This is a list of characters.
    self.alphabet = []
    # emission probability of the I states, one dictionary for each I in the model
    self.eI = [] 
    # emission probability of the M states, one dictionary for each M in the model
    # the first M state is called 'B'; it never emits anything so the associated dictionary is always empty
    self.eM = [{}] 
    # transition probability, one dictionary for each set of states (D,M,I) in the model
    self.t = [] 
    # number of matching states in the model, excluding B and E
    self.nstate = 0
    '''
    def Viterbi(self,query):
        # implement the Viterbi algorithm to find the most probable path of a query sequence and its likelihood
        # return the log-likelihood and the path. If the likelihood is 0, return "$" as the path 
        n = len(query)
        m = self.nstate
        
        VM = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        VI = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        VD = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        
        #outter list length is i/corresponding to the position in the sequence
        
        Vscore = -INF
        aln = "$"
        #VM[i][j]
        # initialization
        VM[0][1] = self.t[0]['MM'] + self.eM[1][query[0]] # B->M1
        VD[0][1] = self.t[0]['MD']                        # B->D1
        VI[0][0] = self.t[0]['MI'] + self.eI[0][query[0]] # B->I0
        
        # YOUR CODE HERE
        #pathVM[0][1].append(0)
        #pathVD[0][1].append(2)
        #pathVI[0][0].append(1)
        
        pathVM = [[list() for j in range(m+2)] for i in range(n+1)]
        pathVI = [[list() for j in range(m+2)] for i in range(n+1)]
        pathVD = [[list() for j in range(m+2)] for i in range(n+1)]
        
        #initialize VD first
        for j in range(2,m+2):
            VD[0][j] = VD[0][j-1]+self.t[j-1]['DD']
            pathVD[0][j].append(2)
        
        #initialize VM
        for j in range(2,m+1):
            VM[0][j] = self.eM[j][query[0]] + VD[0][j-1]+self.t[j-1]['DM']
            pathVM[0][j].append(2)
        
        #intitialize VI
        #first row
        for j in range(1,m+1):
            VI[0][j] = self.eI[j][query[0]] + VD[0][j] + self.t[j]['DI']
            pathVI[0][j].append(2)
            
        for i in range(1,n):
            VI[i][0] = self.eI[0][query[i]] + VI[i-1][0]+ self.t[0]['II']
            pathVI[i][0].append(1)
        
        #i = n
        #VI[i][0] = VI[i-1][0]+ self.t[0]['II']
        #inner list length is j/corresponding to the number of states
        #outter list length is i/corresponding to the position in the sequence
        
        #j is  model state; j is horozontal
        #i is the point in the subsequence that we are at; i is vertical
        a=0
        for j in range(1,m+1):
            for i in range(1,n):
                
                maxVM = [VM[i-1][j-1] + self.t[j-1]['MM'], VI[i-1][j-1] + self.t[j-1]['IM'], VD[i][j-1] + self.t[j-1]['DM'] ]
                VM[i][j] = self.eM[j][query[i]]+ max(maxVM)
                indic = maxVM.index(max(maxVM))
                if indic == 0: 
                    tempVM = pathVM[i-1][j-1][:]
                if indic == 1: 
                    tempVM = pathVI[i-1][j-1][:]
                if indic == 2:
                    tempVM = pathVD[i][j-1][:]
                tempVM.append(indic)
                
                
                maxVI = [VM[i-1][j] + self.t[j]['MI'], VI[i-1][j] + self.t[j]['II'], VD[i][j] + self.t[j]['DI']]
                VI[i][j] = self.eI[j][query[i]]+ max(maxVI)
                indic = maxVI.index(max(maxVI))
                if indic == 0: 
                    tempVI = pathVM[i-1][j][:]
                if indic == 1: 
                    tempVI = pathVI[i-1][j][:]
                if indic == 2:
                    tempVI = pathVD[i][j][:]
                tempVI.append(indic)
                
                maxVD = [VM[i-1][j-1] + self.t[j-1]['MD'], VI[i-1][j-1] + self.t[j-1]['ID'], VD[i][j-1] + self.t[j-1]['DD']]
                VD[i][j] = max(maxVD)
                indic = maxVD.index(max(maxVD))
                if indic == 0: 
                    tempVD = pathVM[i-1][j-1][:]
                if indic == 1: 
                    tempVD = pathVI[i-1][j-1][:]
                if indic == 2:
                    tempVD = pathVD[i][j-1][:]
                tempVD.append(indic)
                
                pathVM[i][j] = tempVM[:]
                pathVI[i][j] = tempVI[:]
                pathVD[i][j] = tempVD[:]
                
                
                
        #last row of VD
        i = n
        for j in range(0, m+1):
            maxVD = [VM[i-1][j-1] + self.t[j-1]['MD'], VI[i-1][j-1] + self.t[j-1]['ID'], VD[i][j-1] + self.t[j-1]['DD']]
            VD[i][j] = max(maxVD)
            indic = maxVD.index(max(maxVD))
            if indic == 0: 
                pathVD[i][j] = pathVM[i-1][j-1][:]
            if indic == 1: 
                pathVD[i][j] = pathVI[i-1][j-1][:]
            if indic == 2:
                pathVD[i][j] = pathVD[i][j-1][:]
            pathVD[i][j].append(indic)
        
        i=n
        j=m+1
        maxVM = [VM[i-1][j-1] + self.t[j-1]['MM'], VI[i-1][j-1] + self.t[j-1]['IM'], VD[i][j-1] + self.t[j-1]['DM'] ]
        VM[i][j] = max(maxVM)
        
        indic = maxVM.index(max(maxVM))
        if indic == 0: 
            pathVM[i][j] = pathVM[i-1][j-1][:]
        if indic == 1: 
            pathVM[i][j] = pathVI[i-1][j-1][:]
        if indic == 2:
            pathVM[i][j] = pathVD[i][j-1][:]
        pathVM[i][j].append(indic)
        
        if VM[i][j] != -INF:
            aln = ""
            l = 0
            for k in pathVM[i][j]:
                if k == 0:
                    aln = aln + (query[l].upper())
                    l = l+1
                elif k == 1:
                    aln = aln + (query[l].lower())
                    l = l+1
                elif k == 2:
                    aln = aln + '-'
        
        Vscore = VM[n][m+1]
        return Vscore,aln        


    def Forward(self,query):
        # implement the Forward algorithm to compute the marginal probability of a query sequence
        n = len(query)
        m = self.nstate
        
        FM = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        FI = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        FD = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        
        # initialization
        FM[0][1] = self.t[0]['MM'] + self.eM[1][query[0]]  # B->M1
        FD[0][1] = self.t[0]['MD']                         # B->D1
        FI[0][0] = self.t[0]['MI'] + self.eI[0][query[0]]  # B->I0

        # YOUR CODE HERE
        #pathFM[0][1].append(0)
        #pathFD[0][1].append(2)
        #pathFI[0][0].append(1)
        '''
        pathFM = [[list() for j in range(m+2)] for i in range(n+1)]
        pathFI = [[list() for j in range(m+2)] for i in range(n+1)]
        pathFD = [[list() for j in range(m+2)] for i in range(n+1)]
        '''
        #initialize FD first
        for j in range(2,m+2):
            FD[0][j] = FD[0][j-1]+self.t[j-1]['DD']
            #pathFD[0][j].append(2)
        
        #initialize FM
        for j in range(2,m+1):
            FM[0][j] = self.eM[j][query[0]] + FD[0][j-1]+self.t[j-1]['DM']
            #pathFM[0][j].append(2)
        
        #intitialize FI
        #first row
        for j in range(1,m+1):
            FI[0][j] = self.eI[j][query[0]] + FD[0][j] + self.t[j]['DI']
            #pathFI[0][j].append(2)
            
        for i in range(1,n):
            FI[i][0] = self.eI[0][query[i]] + FI[i-1][0]+ self.t[0]['II']
            #pathFI[i][0].append(1)
        
        #i = n
        #FI[i][0] = FI[i-1][0]+ self.t[0]['II']
        #inner list length is j/corresponding to the number of states
        #outter list length is i/corresponding to the position in the sequence
        
        #j is  model state; j is horozontal
        #i is the point in the subsequence that we are at; i is vertical
        a=0
        for j in range(1,m+1):
            for i in range(1,n):
                
                maxFM = [FM[i-1][j-1] + self.t[j-1]['MM'], FI[i-1][j-1] + self.t[j-1]['IM'], FD[i][j-1] + self.t[j-1]['DM'] ]
                FM[i][j] = self.eM[j][query[i]]+ log_sum_exp(maxFM)
                '''
                indic = maxFM.index(max(maxFM))
                if indic == 0: 
                    tempFM = pathFM[i-1][j-1][:]
                if indic == 1: 
                    tempFM = pathFI[i-1][j-1][:]
                if indic == 2:
                    tempFM = pathFD[i][j-1][:]
                tempFM.append(indic)
                '''
                
                maxFI = [FM[i-1][j] + self.t[j]['MI'], FI[i-1][j] + self.t[j]['II'], FD[i][j] + self.t[j]['DI']]
                FI[i][j] = self.eI[j][query[i]]+ log_sum_exp(maxFI)
                '''
                indic = maxFI.index(max(maxFI))
                if indic == 0: 
                    tempFI = pathFM[i-1][j][:]
                if indic == 1: 
                    tempFI = pathFI[i-1][j][:]
                if indic == 2:
                    tempFI = pathFD[i][j][:]
                tempFI.append(indic)
                '''
                
                maxFD = [FM[i-1][j-1] + self.t[j-1]['MD'], FI[i-1][j-1] + self.t[j-1]['ID'], FD[i][j-1] + self.t[j-1]['DD']]
                FD[i][j] = log_sum_exp(maxFD)
                '''
                indic = maxFD.index(max(maxFD))
                if indic == 0: 
                    tempFD = pathFM[i-1][j-1][:]
                if indic == 1: 
                    tempFD = pathFI[i-1][j-1][:]
                if indic == 2:
                    tempFD = pathFD[i][j-1][:]
                tempFD.append(indic)
                
                pathFM[i][j] = tempFM[:]
                pathFI[i][j] = tempFI[:]
                pathFD[i][j] = tempFD[:]
                '''
                #last row of FD
        i = n
        for j in range(0, m+1):
            maxFD = [FM[i-1][j-1] + self.t[j-1]['MD'], FI[i-1][j-1] + self.t[j-1]['ID'], FD[i][j-1] + self.t[j-1]['DD']]
            FD[i][j] = log_sum_exp(maxFD)
            '''
            indic = maxFD.index(max(maxFD))
            if indic == 0: 
                pathFD[i][j] = pathFM[i-1][j-1][:]
            if indic == 1: 
                pathFD[i][j] = pathFI[i-1][j-1][:]
            if indic == 2:
                pathFD[i][j] = pathFD[i][j-1][:]
            pathFD[i][j].append(indic)
            '''
        
        i=n
        j=m+1
        maxFM = [FM[i-1][j-1] + self.t[j-1]['MM'], FI[i-1][j-1] + self.t[j-1]['IM'], FD[i][j-1] + self.t[j-1]['DM'] ]
        FM[i][j] = log_sum_exp(maxFM)
        '''
        indic = maxFM.index(max(maxFM))
        if indic == 0: 
            pathFM[i][j] = pathFM[i-1][j-1][:]
        if indic == 1: 
            pathFM[i][j] = pathFI[i-1][j-1][:]
        if indic == 2:
            pathFM[i][j] = pathFD[i][j-1][:]
        pathFM[i][j].append(indic)
        
        if FM[i][j] != -INF:
            aln = ""
            l = 0
            for k in pathFM[i][j]:
                if k == 0:
                    aln = aln + (query[l].upper())
                    l = l+1
                elif k == 1:
                    aln = aln + (query[l].lower())
                    l = l+1
                elif k == 2:
                    aln = aln + '-'
        '''
        Fscore = FM[n][m+1]        
        return Fscore