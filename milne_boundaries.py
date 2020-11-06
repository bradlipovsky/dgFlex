from dolfin import *

def all_boundaries_clamped(x, on_boundary):
    
    # Define several regions of the boundary.
    UpperRight  = (x[0] > -465.7e3) & (x[1] > -630.0e3) & (x[1] < -626.65e3)
    MiddleRight = (x[0] > -466.2e3) & (x[1] > -636.0e3) & (x[1] < -631.0e3)
    #MiddleRight = (x[0] > -466.2e3) & (x[1] > -635.0e3) & (x[1] < -630.9e3)
    LowerRight  = (x[0] > -468.5e3) & (x[1] < -636.0e3)
    UpperLeft   = (x[0] < -475.0e3) & (x[1] < -635.5e3) & (x[1] > -642e3)
    MiddleLeft  = (x[0] < -470.0e3) & (x[1] < -638.5e3)
    
    return on_boundary and (MiddleLeft or UpperRight or UpperLeft or MiddleRight or LowerRight)

def weakened_upper_left(x, on_boundary):
    
    # Define several regions of the boundary.
    UpperRight  = (x[0] > -465.7e3) & (x[1] > -630.0e3) & (x[1] < -626.65e3)
    MiddleRight = (x[0] > -466.2e3) & (x[1] > -636.0e3) & (x[1] < -631.0e3)
    LowerRight  = (x[0] > -468.1e3) & (x[1] < -636.0e3)
    UpperLeft   = (x[0] < -475.0e3) & (x[1] < -635.5e3) & (x[1] > -642e3)
    MiddleLeft  = (x[0] < -473.0e3) & (x[1] < -638.5e3)
    
    return on_boundary and (MiddleLeft or UpperRight or MiddleRight or LowerRight)

def weakened_uppers(x, on_boundary):
    
    # Define several regions of the boundary.
    UpperRight  = (x[0] > -465.7e3) & (x[1] > -630.0e3) & (x[1] < -626.65e3)
    MiddleRight = (x[0] > -466.2e3) & (x[1] > -636.0e3) & (x[1] < -631.0e3)
    LowerRight  = (x[0] > -468.1e3) & (x[1] < -636.0e3)
    UpperLeft   = (x[0] < -475.0e3) & (x[1] < -635.5e3) & (x[1] > -642e3)
    MiddleLeft  = (x[0] < -473.0e3) & (x[1] < -638.5e3)
    
    return on_boundary and (MiddleLeft or MiddleRight or LowerRight)

class upper_left_slip(SubDomain):
	def inside(self, x, on_boundary):
		UpperLeft   = (x[0] < -475.0e3) & (x[1] < -635.5e3) & (x[1] > -642e3)
		return on_boundary and UpperLeft

class upper_left_and_right_slip(SubDomain):
	def inside(self, x, on_boundary):
		UpperRight  = (x[0] > -465.7e3) & (x[1] > -630.0e3) & (x[1] < -626.65e3)
		UpperLeft   = (x[0] < -475.0e3) & (x[1] < -635.5e3) & (x[1] > -642e3)
		return on_boundary and UpperLeft and UpperRight

class clamped_margins(SubDomain):
	def inside(self, x, on_boundary):
    
	    # Define several regions of the boundary.
		UpperRight  = (x[0] > -465.7e3) & (x[1] > -630.0e3) & (x[1] < -626.65e3)
		MiddleRight = (x[0] > -466.2e3) & (x[1] > -636.0e3) & (x[1] < -631.0e3)
	    #MiddleRight = (x[0] > -466.2e3) & (x[1] > -635.0e3) & (x[1] < -630.9e3)
		LowerRight  = (x[0] > -468.5e3) & (x[1] < -636.0e3)
		UpperLeft   = (x[0] < -475.0e3) & (x[1] < -635.5e3) & (x[1] > -642e3)
		MiddleLeft  = (x[0] < -470.0e3) & (x[1] < -638.5e3)
		return on_boundary and (MiddleLeft or UpperRight or UpperLeft or MiddleRight or LowerRight)

