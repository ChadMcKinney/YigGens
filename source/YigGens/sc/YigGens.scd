// YigClifford, an implementation of A Clifford (Pickover) Strange Attractor

/////////////////
// YigCliffordN
/////////////////

// Noisey
{ YigCliffordN.ar(MouseX.kr(20,20000,1), MouseY.kr(-2, 2)).poll;  }.play;

// Sci-fi
{ SinOsc.ar(YigCliffordN.ar(MouseX.kr(1,2000,1), MouseY.kr(-2, 2)).madd(0.5,0.5).poll.linexp(0, 1, 20, 20000), 0, 0.1)  }.play;

// Glitchy
{ y=YigCliffordN.ar(MouseX.kr(1,2000,1)).madd(0.5,0.5);BufCombC.ar(LocalBuf(44100), SinOsc.ar(MouseY.kr(1, 44100, 1),0,0.1),1/y.linexp(0, 1, 1, 20000),5-(y*5)) }.play

// Drumming
{ GVerb.ar(Decay.ar(LeakDC.ar(YigCliffordN.ar(10, MouseX.kr(-2, -0.3, 0), MouseY.kr(2.5, 3, 0), MouseX.kr(0,1), MouseY.kr(0, 1))/3), 3), 10, 1.5, 0.3)  }.play;


/////////////////
// YigCliffordL
/////////////////

// Noisey
{ YigCliffordL.ar(MouseX.kr(20,20000,1), MouseY.kr(-2, 2)).poll;  }.play;

// Sci-fi
{ SinOsc.ar(YigCliffordL.ar(MouseX.kr(1,2000,1), MouseY.kr(-2, 2)).madd(0.5,0.5).poll.linexp(0, 1, 20, 20000), 0, 0.1)  }.play;

// Glitchy
{ y=YigCliffordL.ar(MouseX.kr(1,2000,1), LFNoise1.kr(1,2)).madd(0.5,0.5);BufCombC.ar(LocalBuf(44100), SinOsc.ar(MouseY.kr(1, 44100, 1),0,0.1),1/y.linexp(0, 1, 1, 20000),5-(y*5)) }.play

// Drumming
{ GVerb.ar(Decay.ar(LeakDC.ar(YigCliffordL.ar(10, MouseX.kr(-2, -0.3, 0), MouseY.kr(2.5, 3, 0), MouseX.kr(0,1), MouseY.kr(0, 1))/3), 3), 10, 1.5, 0.3)  }.play;

/////////////////
// YigCliffordC
/////////////////

// Noisey
{ YigCliffordC.ar(MouseX.kr(20,20000,1), MouseY.kr(-2, 2)).poll;  }.play;

// Sci-fi
{ SinOsc.ar(YigCliffordC.ar(MouseX.kr(1,2000,1), MouseY.kr(-2, 2)).madd(0.5,0.5).poll.linexp(0, 1, 20, 20000), 0, 0.1)  }.play;

// Glitchy
{ y=YigCliffordL.ar(MouseX.kr(1,2000,1), LFNoise1.kr(1,2)).madd(0.5,0.5);BufCombC.ar(LocalBuf(44100), SinOsc.ar(MouseY.kr(1, 44100, 1),0,0.1),1/y.linexp(0, 1, 1, 20000),5-(y*5)) }.play


//////////////////
// YigClifford3DN
//////////////////

/* Remember that this outputs on 3 channels, [x, y, z]. If you want to use all three channels you need to explicitly index for the 3rd channel for the output to go to your speakers */

// Noisey
({
	y = YigClifford3DN.ar(MouseX.kr(20,20000,1)).poll;
	Out.ar(0, [y[0]+(y[2]/2), y[1] + (y[2]/2)]/2);
}.play)

// Sci-fi
({
	y = SinOsc.ar(
		YigClifford3DN.ar(
			MouseX.kr(1,2000,1),
			MouseY.kr(-2, 2)
		).madd(0.5,0.5).poll.linexp(0, 1, 20, 20000),
		0,
		0.1);
	Out.ar(0, [y[0]+(y[2]/2), y[1] + (y[2]/2)]);
}.play)

// Glitchy
({
	y = YigClifford3DN.ar(MouseX.kr(1,2000,1)).madd(0.5,0.5);
	b = BufCombC.ar(LocalBuf(44100), SinOsc.ar(MouseY.kr(1, 44100, 1),0,0.1),1/y.linexp(0, 1, 1, 20000),5-(y*5));
	Out.ar(0, [b[0]+(b[2]/2), b[1] + (b[2]/2)]/2);
}.play)

// Drumming
({
	y = YigClifford3DN.ar(10, MouseX.kr(-2, -0.3, 0), MouseY.kr(2.5, 3, 0), MouseX.kr(0,1), MouseY.kr(0, 1));
	GVerb.ar(
		Decay.ar(
			LeakDC.ar([y[0]+(y[2]/2), y[1] + (y[2]/2)]/3),
			3),
		10,
		1.5,
		0.3);

}.play)


//////////////////
// YigClifford3DL
//////////////////

/* Remember that this outputs on 3 channels, [x, y, z]. If you want to use all three channels you need to explicitly index for the 3rd channel for the output to go to your speakers */

// Noisey
({
	y = YigClifford3DL.ar(MouseX.kr(20,20000,1)).poll;
	Out.ar(0, [y[0]+(y[2]/2), y[1] + (y[2]/2)]/2);
}.play)

// Sci-fi
({
	y = SinOsc.ar(
		YigClifford3DL.ar(
			MouseX.kr(1,2000,1),
			MouseY.kr(-2, 2)
		).madd(0.5,0.5).poll.linexp(0, 1, 20, 20000),
		0,
		0.1);
	Out.ar(0, [y[0]+(y[2]/2), y[1] + (y[2]/2)]);
}.play)

// Glitchy
({
	y = YigClifford3DL.ar(MouseX.kr(1,2000,1)).madd(0.5,0.5);
	b = BufCombC.ar(LocalBuf(44100), SinOsc.ar(MouseY.kr(1, 44100, 1),0,0.1),1/y.linexp(0, 1, 1, 20000),5-(y*5));
	Out.ar(0, [b[0]+(b[2]/2), b[1] + (b[2]/2)]/2);
}.play)

// Drumming
({
	y = YigClifford3DL.ar(10, MouseX.kr(-2, -0.3, 0), MouseY.kr(2.5, 3, 0), MouseX.kr(0,1), MouseY.kr(0, 1));
	GVerb.ar(
		Decay.ar(
			LeakDC.ar([y[0]+(y[2]/2), y[1] + (y[2]/2)]/3),
			3),
		10,
		1.5,
		0.3);

}.play)


//////////////////
// YigClifford3DC
//////////////////

/* Remember that this outputs on 3 channels, [x, y, z]. If you want to use all three channels you need to explicitly index for the 3rd channel for the output to go to your speakers */

// Noisey
({
	y = YigClifford3DC.ar(MouseX.kr(20,20000,1)).poll;
	Out.ar(0, [y[0]+(y[2]/2), y[1] + (y[2]/2)]/2);
}.play)

// Sci-fi
({
	y = SinOsc.ar(
		YigClifford3DC.ar(
			MouseX.kr(1,2000,1),
			MouseY.kr(-2, 2)
		).madd(0.5,0.5).poll.linexp(0, 1, 20, 20000),
		0,
		0.1);
	Out.ar(0, [y[0]+(y[2]/2), y[1] + (y[2]/2)]);
}.play)

// Glitchy
({
	y = YigClifford3DC.ar(MouseX.kr(1,2000,1)).madd(0.5,0.5);
	b = BufCombC.ar(LocalBuf(44100), SinOsc.ar(MouseY.kr(1, 44100, 1),0,0.1),1/y.linexp(0, 1, 1, 20000),5-(y*5));
	Out.ar(0, [b[0]+(b[2]/2), b[1] + (b[2]/2)]/2);
}.play)

//////////////////
// YigMandelbulbN
//////////////////

// Noisey
({
	y = LeakDC.ar(YigMandelbulbN.ar(MouseY.kr(1, 2000, 1), 8, 8, MouseX.kr(1, 1,32768).round)).softclip.poll;
	Out.ar(0, [y[0]+(y[2]/2), y[1] + (y[2]/2)]/2);
}.play )

// Sci-fi
({
	m = YigMandelbulbN.ar(
			MouseY.kr(1, 1000, 1),
			MouseX.kr(3, 12).round,
			MouseY.kr(3, 8),
			MouseX.kr(1, 32768).round
		).madd(0.5,0.5).poll;
	y = SinOsc.ar(
		m[0].linexp(0, 1, 1, 10000) + (SinOsc.ar(m[1].linexp(0, 1, 1, 200),0) * m[2].linexp(0, 1, 1, 10000)),
		SinOsc.ar(m * 10 + 10, 0, 0.5, 0.5),
		0.125);
	Out.ar(0, [y[0]+(y[2]/2), y[1] + (y[2]/2)]);
}.play)

// Actually Glitchy
({
	y = YigMandelbulbN.ar(MouseX.kr(20,4000,1), MouseX.kr(3, 16,1), MouseY.kr(4, 32).round, MouseX.kr(1,32768).round);
	b = BufCombC.ar(LocalBuf(44100), SinOsc.ar(MouseY.kr(1, 44100, 1),0,0.1),1/y.softclip.linexp(0, 1, 1, 20000),5-(y*5));
	Out.ar(0, [b[0]+(b[2]/2), b[1] + (b[2]/2)]/2);
}.play)

// Drumming
({
	y = LeakDC.ar(YigMandelbulbN.ar(200, MouseY.kr(3, 11).round, MouseY.kr(4, 10).round, MouseX.kr(1, 32768).round).softclip);
	y = [BufDelayC.ar(LocalBuf(44100 * 4), y[0] + (y[2]/2), y[1].linlin(-1, 1, 1, 6).round, 1, y),
		 BufDelayC.ar(LocalBuf(44100 * 4), y[1] + (y[2]/2), y[0].linlin(-1, 1, 1, 6).round, 1, y)];
	LPF.ar(GVerb.ar(
		BLowShelf.ar(LeakDC.ar(y) * 5, 80, 12),
		30,
		2,
		0.1,
		0.1), 3000) * 0.1;

}.play