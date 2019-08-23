class algorithm(object):
    """
    """
    def __init__(self,Tag,ParameterStartRange,ParameterEndRange,ParameterScaleFactor,ParameterScaleOperator,IndirectIndexFunc,SpecialFunc,ScaleFactorList,ScaleOperatorList,NodeStartOffsetList=[0]):
        """
	"""
	self.Tag=Tag
	self.NumParameters=len(ParameterStartRange)
	self.InputParameterStartRange=ParameterStartRange
	self.InputParameterEndRange=ParameterEndRange
	self.InputParameterScaleFactor=ParameterScaleFactor
	self.InputParameterScaleOperator=ParameterScaleOperator
        self.ScaleFactorList=ScaleFactorList
	self.ScaleOperatorList=ScaleOperatorList
        self.IndirectIndexFunc = IndirectIndexFunc
        self.SpecialFunc = SpecialFunc
        self.NodeStartOffsetList=NodeStartList
        self.TagList=[]
	for i in range(len( self.InputParameterStartRange)):
	    if (self.InputParameterStartRange[i] < self.InputParameterEndRange[i]):
	        self.TagList.append(-1)
	    elif (self.InputParameterStartRange[i] > self.InputParameterEndRange[i]):
	        self.TagList.append(1)
	    else:
	        self.TagList.append(0)

    def scale(self,Parameters,Index):
        """
        """
	for i in range(len(Parameters)):
	    Parameters[i] = self.ScaleOperatorList[self.IndirectIndexFunc(Index)][i](Parameters[i],self.ScaleFactorList[self.IndirectIndexFunc(Index)][i])
