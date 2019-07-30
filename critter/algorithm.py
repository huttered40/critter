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
        self.Save = 0

    def next(self):
        """
        """
	if (self.Save == 1):
	    return 0
	elif (self.CurrentStartParameters == self.InputParameterEndRange):
	    self.Save=1
	    return 1
        else:
	    for i in range(len(self.InputParameterStartRange)):
                self.CurrentStartParameters[i] = self.InputParameterScaleOperator[i](self.CurrentStartParameters[i],self.InputParameterScaleFactor[i])
                self.CurrentScaleParameters[i] = self.CurrentStartParameters[i]
	    return 1

    def scale(self,index):
        """
        """
	for i in range(len(self.InputParameterStartRange)):
	    self.CurrentScaleParameters[i] = self.ScaleOperatorList[index][i](self.CurrentScaleParameters[i],self.ScaleFactorList[index][i])
