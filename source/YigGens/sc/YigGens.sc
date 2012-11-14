YigGen : MultiOutUGen {
	
}

YigCliffordN : YigGen {
	*ar { arg freq = 440.0, a = -0.966918, b = 2.879879, c = 0.765145, 
		d = 0.744728, x = 0.01, y = 0.01, mul = 1.0, add = 0.0;
		^this.multiNew('audio', freq, a, b, c, d, x, y).madd(mul, add);
	}
	
	init { arg ... theInputs;
		inputs = theInputs;
		channels = [ 
			OutputProxy(rate, this, 0), 
			OutputProxy(rate, this, 1)
		];
		^channels;
	}
}

YigCliffordL : YigCliffordN {}

YigCliffordC : YigCliffordL {}

Yig3DGen : MultiOutUGen {
	
}

YigClifford3DN : Yig3DGen {
	*ar { arg freq = 440.0, a = 2.24, b = 0.43, c = -0.65, 
		d = -2.43, x = 0.01, y = 0.01, z = 0.01, mul = 1.0, add = 0.0;
		^this.multiNew('audio', freq, a, b, c, d, x, y, z).madd(mul, add);
	}
	
	init { arg ... theInputs;
		inputs = theInputs;
		channels = [ 
			OutputProxy(rate, this, 0), 
			OutputProxy(rate, this, 1),
			OutputProxy(rate, this, 2)
		];
		^channels;
	}
}

YigClifford3DL : YigClifford3DN {}

YigClifford3DC : YigClifford3DL {}

YigMandelbulbN : Yig3DGen {
	// rule is 1 - 32768
	*ar { arg freq = 440.0, p = 8, size = 8, rule = 50, maxIter = 2056, mul = 1.0, add = 0.0;
		^this.multiNew('audio', freq, p, size, rule, maxIter).madd(mul, add);
	}
	
	init { arg ... theInputs;
		inputs = theInputs;
		channels = [ 
			OutputProxy(rate, this, 0), 
			OutputProxy(rate, this, 1),
			OutputProxy(rate, this, 2)
		];
		^channels;
	}
}
