class algorithm(object):
    """
    """
    def __init__(self,Tag,ParameterStartRange,ParameterEndRange,ParameterScaleFactor,ParameterScaleOperator,ScaleFactorList,ScaleOperatorList):
        """
	"""
	self.Tag=Tag
	self.InputParameterStartRange=ParameterStartRange
	self.InputParameterEndRange=ParameterEndRange
	self.InputParameterScaleFactor=ParameterScaleFactor
	self.InputParameterScaleOperator=ParameterScaleOperator
        self.ScaleFactorList=ScaleFactorList
	self.ScaleOperatorList=ScaleOperatorList
	self.CurrentStartParameters=list(self.InputParameterStartRange)
	self.CurrentScaleParameters=list(self.InputParameterStartRange)

    def getNextStartParameters(self):
        """
        """
        for i in range(len(self.InputParameterStartRange)):
            self.CurrentStartParameters[i] = self.InputParameterScaleOperator[i](self.CurrentStartParameters[i],self.InputParameterScaleFactor[i])
            self.CurrentScaleParameters[i] = self.CurrentStartParameters[i]

    def scale(self,index):
        """
        """
	for i in range(len(self.ParameterStartRange)):
	    self.CurrentScaleParameters[i] = self.ScaleOperator[index][i](self.CurrentScaleParameters[i],self.ScaleFactor[index][i])