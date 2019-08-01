class algorithm(object):
    """
    """
    def __init__(self,Tag,ParameterStartRange,ParameterEndRange,ParameterScaleFactor,ParameterScaleOperator,IndirectIndexFunc,SpecialFunc,ScaleFactorList,ScaleOperatorList):
        """
	"""
	self.Tag=Tag
	self.InputParameterStartRange=ParameterStartRange
	self.InputParameterEndRange=ParameterEndRange
	self.InputParameterScaleFactor=ParameterScaleFactor
	self.InputParameterScaleOperator=ParameterScaleOperator
        self.ScaleFactorList=ScaleFactorList
	self.ScaleOperatorList=ScaleOperatorList
        self.IndirectIndexFunc = IndirectIndexFunc
        self.SpecialFunc = SpecialFunc
	self.CurrentStartParameters=list(self.InputParameterStartRange)
	self.CurrentScaleParameters=list(self.InputParameterStartRange)
        self.SaveParameterStartRange=list(self.InputParameterStartRange)
        self.ParameterIndex=0

    def next(self):
        """
        """
        #print("Here in next with ParameterIndex - %d, with val - %d and ref - %d\n"%(self.ParameterIndex, self.CurrentStartParameters[self.ParameterIndex], self.SaveParameterStartRange[self.ParameterIndex]))
        self.CurrentScaleParameters = list(self.CurrentStartParameters)
        while (self.ParameterIndex < len(self.InputParameterStartRange)):
	    while (self.CurrentStartParameters[self.ParameterIndex] != self.InputParameterEndRange[self.ParameterIndex]):
		self.CurrentStartParameters[self.ParameterIndex] = self.InputParameterScaleOperator[self.ParameterIndex](\
		       self.CurrentStartParameters[self.ParameterIndex],self.InputParameterScaleFactor[self.ParameterIndex])
		self.CurrentScaleParameters[self.ParameterIndex] = self.CurrentStartParameters[self.ParameterIndex]
		return 1
	    self.ParameterIndex = self.ParameterIndex + 1
	return 0

    def scale(self,index):
        """
        """
	for i in range(len(self.InputParameterStartRange)):
	    self.CurrentScaleParameters[i] = self.ScaleOperatorList[self.IndirectIndexFunc(index)][i](self.CurrentScaleParameters[i],self.ScaleFactorList[self.IndirectIndexFunc(index)][i])
