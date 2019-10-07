import time
import visa
import numpy
import warnings
import matplotlib
import matplotlib.pyplot as plt
import threading
import PyQt4
import scipy
import scipy.optimize
from PyQt4 import QtGui, QtCore
import sys

class EG_G_7265_DSP_LockIn(object):
    def __init__(self, GPIB_Address = 12, RemoteOnly = False):
        self.VI = visa.instrument('GPIB0::' + str(GPIB_Address))
        self.RemoteOnly(RemoteOnly)
    def __del__(self):
        self.VI.close()
        
    def __chkFloat(self, s):
        if s[-1] == '\x00':
            s = s[:-1]
        return float(s)
        
    def RemoteOnly(self, rO = True):
        if rO:
            self.VI.write('REMOTE 1')
        else:
            self.VI.write('REMOTE 0')
            
    def setTC(self, TC):
        """ Sets the Filter Time Constant
        Usage :
            setTC(Value or Code)
                Codes :
                 '0'  = 10 us
                 '1'  = 20 us
                 '2'  = 40 us
                 '3'  = 80 us
                 '4'  = 160 us
                 '5'  = 320 us
                 '6'  = 640 us
                 '7'  = 5 ms
                 '8'  = 10 ms
                 '9'  = 20 ms
                 '10' = 50 ms
                 '11' = 100 ms
                 '12' = 200 ms
                 '13' = 500 ms
                 '14' = 1 s
                 '15' = 2 s
                 '16' = 5 s
                 '17' = 10 s
                 '18' = 20 s
                 '19' = 50 s
                 '20' = 100 s
                 '21' = 200 s
                 '22' = 500 s
                 '23' = 1 ks
                 '24' = 2 ks
                 '25' = 5 ks
                 '26' = 10 ks
                 '27' = 20 ks
                 '28' = 50 ks
                 '29' = 100 ks
                 '30' = 200 ks
                Values are rounded to corresponding code.
        """
        if type(TC) != type(' '):
            TC = numpy.abs(TC) - 1E-9
            TC = numpy.clip(TC, 0, 199E3)
            bins = numpy.array([10E-6, 20E-6, 40E-6, 80E-6,
                                160E-6, 320E-6, 640E-6,
                                5E-3, 10E-3, 20E-3, 50E-3,
                                100E-3, 200E-3, 500E-3,
                                1, 2, 5, 10, 20, 50,
                                100, 200, 500,
                                1E3, 2E3, 5E3, 10E3, 20E3, 50E3,
                                100E3, 200E3])
            TC = str(len(bins[bins <= TC]))
        if TC in map(str, range(31)):
            self.VI.write('TC ' + TC)
        else:
            warning.warn('EG&G 7265 Lock-In Wrong Time Constant Code')
    def __getTC(self):
        return float(self.VI.ask('TC.'))
    TC = property(__getTC, setTC, None, "Filter Time Constant.")
    
    def setSEN(self, vSen):
        """ Sets the Full Scale Sensitivity
        Usage :
            setTC(Value or Code)
                Codes : IMODE=0 IMODE=1 IMODE=2
                 '1'  = 2 nV    2 fA    n/a
                 '2'  = 5 nV    5 fA    n/a
                 '3'  = 10 nV   10 fA   n/a
                 '4'  = 20 nV   20 fA   n/a
                 '5'  = 50 nV   50 fA   n/a
                 '6'  = 100 nV  100 fA  n/a
                 '7'  = 200 nV  200 fA  2 fA
                 '8'  = 500 nV  500 fA  5 fA
                 '9'  = 1 uV    1 pA    10 fA
                 '10' = 2 uV    2 pA    20 fA
                 '11' = 5 uV    5 pA    50 fA
                 '12' = 10 uV   10 pA   100 fA
                 '13' = 20 uV   20 pA   200 fA
                 '14' = 50 uV   50 pA   500 fA
                 '15' = 100 uV  100 pA  1 pA
                 '16' = 200 uV  200 pA  2 pA
                 '17' = 500 uV  500 pA  5 pA
                 '18' = 1 mV    1 nA    10 pA
                 '19' = 2 mV    2 nA    20 pA
                 '20' = 5 mV    5 nA    50 pA
                 '21' = 10 mV   10 nA   100 pA
                 '22' = 20 mV   20 nA   200 pA
                 '23' = 50 mV   50 nA   500 pA
                 '24' = 100 mV  100 nA  1 nA
                 '25' = 200 mV  200 nA  2 nA
                 '26' = 500 mV  500 nA  5 nA
                 '27' = 1V      1 uA    10 nA
                Values are rounded to corresponding code.
        """
        InputMode = self.VI.ask('IMODE')
        if type(vSen) != type(' '):
            if InputMode == '0':
                vSen = vSen * 1.0E-6 
            vSen = numpy.abs(vSen) - 1E-18
            vSen = numpy.clip(vSen, 0, 0.99E-6)
            bins = numpy.array([0, 2E-15, 5E-15, 10E-15, 20E-15, 
                                50E-15, 100E-15, 200E-15, 500E-15,
                                1E-12, 2E-12, 5E-12, 10E-12, 20E-12, 
                                50E-12, 100E-12, 200E-12, 500E-12,
                                1E-9, 2E-9, 5E-9, 10E-9, 20E-9, 
                                50E-9, 100E-9, 200E-9, 500E-9,
                                1E-6])
            vSen = len(bins[bins <= vSen])
            if InputMode == '2':
                vSen = vSen + 6
                vSen = numpy.clip(vSen, 7, 27)
            vSen = str(vSen)
        if vSen in map(str, range(1,28)):
            self.VI.write('SEN ' + vSen)
        else:
            warning.warn('EG&G 7265 Lock-In Wrong Scale Sensitivity Code')
    def __getSEN(self):
        return float(self.VI.ask('SEN.'))
    SEN = property(__getSEN, setSEN, None, "Full Scale Sensitivity.")
    
    def FilterSlope(self, sl):
        """Set the output filter slope
        Usage :
            FilterSlope(Code)
                Codes :
                 '0' = 6 dB/octave
                 '1' = 12 dB/octave
                 '2' = 18 dB/octave
                 '3' = 24 dB/octave
        """
        if sl in ['0', '1', '2', '3']:
            self.VI.write('SLOPE ' + sl)
        else:
            warning.warn('EG&G 7265 Lock-In Wrong Slope Code')
    
    def InputMode(self, imode):
        """Current/Voltage mode Input Selector
        Usage :
            InputMode(Code)
                Codes :
                 '0' = Voltage Mode
                 '1' = Current Mode High bandwidth
                 '2' = Current Mode Low noise
        """
        if imode in ['0', '1', '2']:
            self.VI.write('IMODE ' + imode)
        else:
            warning.warn('EG&G 7265 Lock-In Wrong Input Mode Code')
            
    def VoltageInputMode(self, vmode):
        """Voltage Input Configuration
        Usage :
            VoltageInputMode(Code)
                Codes :
                 '0' = Grounded
                 '1' = A Input only
                 '2' = -B Input only
                 '3' = A-B diferential mode
        """
        self.VI.write('VMODE 0')
        if vmode in ['0', '1', '2', '3']:
            self.VI.write('VMODE ' + vmode)
        else:
            warning.warn('EG&G 7265 Lock-In Wrong Voltage Configuration Code')
    
    def Sync(self, Sy = True):
        """Enable or disable Synchonous time constant
        """
        if Sy:
            self.VI.write('SYNC 1')
        else:
            self.VI.write('SYNC 0')
            
    def setOscilatorFreq(self, freq):
        """Set the internal Oscilator Frequency"""
        fq = str(int(freq*1E3))
        self.VI.write('OF '+ fq)
        
    def setOscilatorAmp(self, amp):
        """Set the internal Oscilator Amplitude"""
        A = str(int(amp*1E6))
        self.VI.write('OA '+ A)
        
    def setRefPhase(self, ph):
        """Set the phase reference"""
        P = str(int(ph*1E3))
        self.VI.write('REFP '+ P)
        
    def getRefPhase(self):
        """Get the programed phase reference"""
        return self.__chkFloat(self.VI.ask('REFP.'))
        
    def ConfigureInput(self, InDev = 'FET', Coupling = 'AC', Ground = 'GND', AcGain = 'Auto'):
        if InDev  == 'FET':
            self.VI.write('FET 1')
        elif InDev  == 'Bipolar':
            self.VI.write('FET 0')
        else:
            warning.warn('EG&G 7265 Lock-In Wrong Input device control Code')

        if Coupling  == 'DC':
            self.VI.write('CP 1')
        elif Coupling  == 'AC':
            self.VI.write('CP 0')
        else:
            warning.warn('EG&G 7265 Lock-In Wrong Input Coupling Code')

        if Ground  == 'GND':
            self.VI.write('FLOAT 0')
        elif Ground  == 'Float':
            self.VI.write('FLOAT 1')
        else:
            warning.warn('EG&G 7265 Lock-In Wrong Input Coupling Code')
            
        if AcGain == 'Auto':
            self.VI.write('AUTOMATIC 1')
        elif AcGain in map(str, range(10)):
            self.VI.write('AUTOMATIC 0')
            self.VI.write('ACGAIN ' + AcGain)
        else:
            warning.warn('EG&G 7265 Lock-In Wrong AcGain Code')
            
    @property
    def X(self):
        return self.__chkFloat(self.VI.ask('X.'))
    @property
    def Y(self):
        return self.__chkFloat(self.VI.ask('Y.'))
    @property
    def Magnitude(self):
        return self.__chkFloat(self.VI.ask('MAG.'))
    @property
    def Phase(self):
        return self.__chkFloat(self.VI.ask('PHA.'))
    @property
    def Freq(self):
        return self.__chkFloat(self.VI.ask('FRQ.'))
        
class KEPCO_BOP_50_20MG(object):
    def __init__(self, GPIB_Address = 7):
        self.VI = visa.instrument('GPIB0::' + str(GPIB_Address))
        self.VI.write('*CLS')
        self.VI.write('*RST')
        self.VI.write('OUTPUT ON')
    def __del__(self):
        #self.VI.write('CURR 0')
        self.VI.write('VOLT 0')
        self.VI.close()
        
    def __SetRange(self, r):
        """Change operating range for output current or voltage
        Usage :
            Range('Full' / '1/4' / 'AUTO')
        """
        ValidCodes = ['Full' , '1/4' , 'AUTO']
        if r in ValidCodes:
            OutStrings = [
                          ['VOLT:RANG 1', 'CURR:RANG 1'],
                          ['VOLT:RANG 4', 'CURR:RANG 4'],
                          ['VOLT:RANG:AUTO', 'CURR:RANG:AUTO']
                         ]
            self.VI.write(OutStrings[ValidCodes == r][int(self.VI.ask('FUNC:MODE?'))])
        else:
            warnings.warn('Range error code')
        
    def Output(self, out):
        """Enable or disable power supply output
        Usage :
            Output('ON'/'OFF')
        """
        if out in ['ON', 'OFF']:
            self.VI.write('OUTPUT ' + out)
        else:
            warnings.warn('Output error code')
        
    def CurrentMode(self):
        """ Changes to constant current operation mode """
        self.VI.ask('SYST:ERR?')
        self.VI.write('FUNC:MODE CURR')
    
    def VoltageMode(self):
        """ Changes to constant voltage operation mode """
        self.VI.ask('SYST:ERR?')
        self.VI.write('FUNC:MODE VOLT')
    
    @property
    def OperationMode(self):
        """Returns actual operation mode"""
        modes = ['Constant Voltage', 'Constant Current']
        return modes[int(self.VI.ask('FUNC:MODE?'))]
    
    def VoltageOut(self, vOut):
        """ Sets the Output Voltage
        Usage :
            VoltageOut(Voltage)
        """
        self.VI.ask('SYST:ERR?')
        self.VI.write('VOLT ' + "%0.4f" % vOut)
    def __getVoltageOut(self):
        self.VI.ask('SYST:ERR?')
        return float(self.VI.ask('VOLT?'))
    Voltage = property(__getVoltageOut, VoltageOut, None, 
            """ On Voltage mode:
            Sets output voltage and return programed voltage
            On Current mode:
            Sets and return protection voltage
            """)
    
    def CurrentOut(self, cOut):
        """ Sets the Output Current
        Usage :
            CurrentOut(Current)
        """
        self.VI.ask('SYST:ERR?')
        self.VI.write('CURR ' + "%0.4f" % cOut)
    def __getCurrentOut(self):
        self.VI.ask('SYST:ERR?')
        return float(self.VI.ask('CURR?'))
    Current = property(__getCurrentOut, CurrentOut, None,
            """ On Voltage mode:
            Sets and return protection current
            On Current mode:
            Sets output current and return programed current
            """)

    def BEEP(self):
        """BEEP"""
        self.VI.ask('SYST:ERR?')
        self.VI.write('SYST:BEEP')
        time.sleep(1)

    @property
    def MeasuredVoltage(self):
        """Measured Voltage Value"""
        self.VI.ask('SYST:ERR?')
        return float(self.VI.ask('MEAS:VOLT?'))
    @property
    def MeasuredCurrent(self):
        """Measured Current Value"""
        self.VI.ask('SYST:ERR?')
        return float(self.VI.ask('MEAS:CURR?'))

class KEITHLEY_2000_Multimeter(object):
    def __init__(self, GPIB_Address = 16):
        self.VI = visa.instrument('GPIB0::' + str(GPIB_Address), term_chars = visa.LF)
        time.sleep(0.5)
        self.VI.write('*CLS')
        self.VI.write('*RST')
        time.sleep(0.5)
        self.VI.write('FUNC "VOLT:DC"')
        
        self.VI.write('VOLT:DC:RANG:UPP 1')
        self.VI.write('VOLT:DC:AVER:TCON MOV')
        self.VI.write('VOLT:DC:AVER:COUN 20')
        self.VI.write('VOLT:DC:DIG 5')
        self.VI.write('INIT:CONT 1')
        
    @property
    def Voltage(self):
        """Voltage Value"""
        return float(self.VI.ask('FETCH?'))

class FieldControler(object):
    """
    Controlador de Campo Magnetico
    """
    def __init__(self, parent =  None):
        self.Kepco = KEPCO_BOP_50_20MG()
        self.NanoVolt = KEITHLEY_2000_Multimeter()
        
        self.HperIOut = 223.4 #Oe por I de la salida (aprox)
        self.MaxHRate = 100.0 #Rate H maximo en Oe/s
        self.VInToH = 10.0E3 #Oe por Voltios medidos
        
        self.Kepco.CurrentMode()
        self.Kepco.Voltage = 25.0
        
        #self.NanoVolt.Range(2.0)
        
        try:
            self.log = parent.log
        except:
            self.log = _empty_log
        
    def getField(self, Res = 'Fast', Unit = 'Oe'):
        """Returns the measured magnetic field in Oe.
        Usage :
            getField()
            getField(Res , Unit)
            
        Resolution codes:
            'Fast' (Default) : Return the value of 1 measurement.
            'High' : Return the mean of 10 measurements.
        Units:
            'Oe' : Field in Oe
            'A/m' : Field in A/m
        """
        if Res == 'Fast':
            vIn = self.NanoVolt.Voltage
        elif Res == 'High':
            n = 10
            vsIn = numpy.zeros(n)
            for i in range(n):
                vsIn[i] = self.NanoVolt.Voltage
            vIn = vsIn.mean()
        
        if Unit == 'Oe':
            return self.VInToH * vIn
        elif Unit == 'A/m':
            return self.VInToH * vIn * 1.0E3 / (4 * numpy.pi) 

    def test(self, i):
        self.Kepco.Current = i
        c0 = time.clock()
        n = 25
        vsIn = numpy.zeros(n)
        for a in range(n):
            vsIn[a] = self.NanoVolt.Voltage
        c1 = time.clock()
        ts = numpy.linspace(0, c1-c0, n)
        return (vsIn, ts)
            
    def setField(self, Fld, Alg = 'Fast'):
        """
        Set Magnetic Field
        """
        
        self.log('Setting field : %.1f Oe ... ' % Fld, EOL = '')
        if Alg == 'Fast':
            vInTol = 0.5 / self.VInToH
        elif Alg == 'Exact':
            vInTol = 0.5 / self.VInToH
        else:
            vInTol = 5.0 / self.VInToH
        
        trgVin = Fld / self.VInToH
        actVin = self.NanoVolt.Voltage
        vInErr = numpy.abs(trgVin - actVin)

        while vInErr >= vInTol:
            dIOut = 0.90 * (trgVin - actVin) * self.VInToH / self.HperIOut  #voltaje de debe variar la salida
            #Rampa variando dVolt
            iIni = self.Kepco.Current
            iFin = iIni + dIOut
            step = 0.08 #Amps
            dt = step * self.HperIOut / self.MaxHRate 
            iPoints = numpy.arange(iIni, iFin, step * numpy.sign(dIOut))
            if len(iPoints) < 5:
                iPoints = numpy.linspace(iIni, iFin, 5)
            for i in iPoints:
                self.Kepco.Current = i
                #time.sleep(dt) #hasta quitar la comprovacion de errores de la kepco.

            #Espera a que se estabilice la corriente
            time.sleep(0.1)

            #Lee la proyeccion del campo
            ff = lambda p, t: p[0] * (1 - exp(-p[1] * (t + p[2]))) + p[3]
            errFunc = lambda p, t, vi: (ff(p,t) - vi)**2
            

            optOk = 5
            while optOk >= 5:
                c0 = time.clock()
                n = 20
                vsIn = numpy.zeros(n)
                for a in range(n):
                    vsIn[a] = self.NanoVolt.Voltage
                    time.sleep(0.2)
                c1 = time.clock()
                ts = numpy.linspace(0, c1-c0, n)
                p0 = [trgVin - vsIn[0], 2, 0.1, vsIn[0]]
                px = scipy.optimize.leastsq(errFunc, p0, args = (ts,vsIn), warning = False)
                optOk = px[1]
            actVin = px[0][0] + px[0][3]
            
            
            #calula de nuevo el error 
            if numpy.sign(trgVin - actVin) == numpy.sign(dIOut):
                vInErr = numpy.abs(trgVin - actVin)
            else:
                vInErr = 0
        time.sleep(2)
        self.log('Done.', [125,125,125])
    
    def TurnOff(self):
        self.log('Turning field off ... ', EOL = '')
        iIni = self.Kepco.Current
        tRamp = numpy.abs(iIni) / (self.MaxHRate / self.HperIOut)
        iPoints = numpy.linspace(iIni, 0, 200)
        dt = numpy.abs(iIni) * self.HperIOut  / (self.MaxHRate * 200.0)
        for i in iPoints:
            self.Kepco.Current = i
            #time.sleep(dt) #hasta quitar la comprovacion de errores de la kepco.
        self.log('Done.', [125,125,125])

    def BEEP(self):
        self.Kepco.BEEP()
        

class VSMControler(object):
    """
    Controlador y sensor del VSM
    """
    def __init__(self, parent = None):
        self.LockIn = EG_G_7265_DSP_LockIn(RemoteOnly = False)
        self.LockIn.InputMode('0')
        self.LockIn.VoltageInputMode('1')
        self.LockIn.FilterSlope('3')
        self.LockIn.setRefPhase(-110.635)
        self.confDriver()
        self.confInput()
        #self.emu_per_V = 351.86
        self.emu_per_V = 46.6745
        
        try:
            self.log = parent.log
        except:
            self.log = _empty_log
        
        
    def confDriver(self, OscFrec = 82, OscAmp = 0.63):
        self.LockIn.setOscilatorAmp(OscAmp)
        self.LockIn.setOscilatorFreq(OscFrec)
    def confInput(self, Sen = 1, TC = 0.1, AcGain = '0'):
        self.LockIn.TC = TC
        self.LockIn.SEN = Sen
        self.LockIn.ConfigureInput(AcGain = AcGain)
        
    def ZeroPhase(self):
        TCtemp = self.LockIn.TC
        self.LockIn.TC = 1
        time.sleep(15)
        ph = 0
        for i in range(10):
            time.sleep(1)
            ph = self.LockIn.Phase + ph
        ph = ph / 10.0
        self.LockIn.setRefPhase(self.LockIn.getRefPhase() + ph)
        self.LockIn.TC = TCtemp
        time.sleep(3)
        
    def getRefPhase(self):
        return self.LockIn.getRefPhase()
        
    def getMagnetization(self, n = 20, iniDelay = 1, measDelay = 0, stat = False, tol = 0.05, maxIts = 50):
        self.log('Measuring Magnetization ... ', EOL = '')
        vsIn = numpy.zeros(n)
        time.sleep(iniDelay)
        for i in range(n):
            time.sleep(measDelay)
            vsIn[i] = self.LockIn.X
        vIn = vsIn.mean()
        sigma = vsIn.std()
        maxSigma = numpy.abs(self.LockIn.SEN * tol)

        if stat:
            its = 0
            while (sigma > maxSigma) and (its < maxIts):
                its = its + 1 
                err = (vsIn - vIn)**2
                vsIn = vsIn[err < sigma**2]
                while len(vsIn) < n:
                    time.sleep(measDelay)
                    vsIn = numpy.append(vsIn, self.LockIn.X)
                vIn = vsIn.mean()
                sigma = vsIn.std()
        self.log('Done.', [125,125,125])
        self.log('M = %.3E     | ' % (vIn * self.emu_per_V), [100,100,100], EOL = '')
        self.log('s = %.3E ' % (sigma * self.emu_per_V), [190,190,190])
        return numpy.array([vIn, sigma])* self.emu_per_V
        
    def getAmplitude(self, n = 20, iniDelay = 1, measDelay = 0):
        vsIn = numpy.zeros(n)
        time.sleep(iniDelay)
        for i in range(n):
            time.sleep(measDelay)
            vsIn[i] = self.LockIn.Magnitude
        vIn = vsIn.mean()
        return vIn

class VSM(object):
    class Logger(QtCore.QObject):
        logSignal = QtCore.pyqtSignal(str, QtGui.QColor, int, str)
        def __init__(self):
            QtCore.QObject.__init__(self)
            self.app = QtGui.QApplication(sys.argv)
            self.LogWindow = QtGui.QTextEdit()
            self.LogWindow.setReadOnly(True)
            self.LogWindow.setWindowTitle('VSM Log')
            self.LogWindow.resize(400,650)
            self.LogWindow.setWordWrapMode(0) #No Word Wrap
            self.logSignal.connect(self.logSlot)
        def logSlot(self, log, color, lev, EOL):
            self.LogWindow.moveCursor(QtGui.QTextCursor.End, QtGui.QTextCursor.MoveAnchor) 
            self.LogWindow.setTextColor(color)
            self.LogWindow.insertPlainText(log + EOL)
            self.LogWindow.ensureCursorVisible()
            self.LogWindow.show()
        def log(self, log, color = 'k', level = 0, EOL = '\n'):
            colorDict = {'k':[0,0,0], 'r':[255,0,0], 'g':[0,255,0], 'b':[0,0,255]}
            if color in colorDict.keys():
                color = QtGui.QColor(*colorDict[color])
            else:
                try:
                    color = QtGui.QColor(*color)
                except:
                    color = QtGui.QColor(0,0,0)
            self.logSignal.emit(log, color, level, EOL)

    class MThread(threading.Thread):
        def __init__(self):
            threading.Thread.__init__(self)
            self.Stop = False
            self.clear()
        def __del__(self):
            self.Stop = True
        def addFunction(self, funct, ret = None, *args, **kwargs):
            self.functionStack.append(funct)
            self.argsStack.append(args)
            self.kwargsStack.append(kwargs)
            self.returnStack.append(ret)
        def clear(self):
            self.argsStack = []
            self.kwargsStack = []
            self.functionStack = []
            self.returnStack = []
            self.Out = None
        def run(self):
            while not(self.Stop):
                if len(self.functionStack) == 0:
                    PyQt4.QtTest.QTest.qWait(50)
                else:
                    self.Out = self.functionStack.pop(0)(*self.argsStack.pop(0), **self.kwargsStack.pop(0))
                    try:
                        retF = self.returnStack.pop(0)
                        if callable(retF):
                            retF(self.Out)
                    except:
                        pass
                    
    class HData():
        def __init__(self):
            self.reset()
        def addPoint(self, C = numpy.NaN, M  = numpy.NaN, S  = numpy.NaN):
            pt = numpy.array([[C],[M],[S]])
            isnanpt = numpy.isnan(pt)[:,-1]
            isnandat = numpy.isnan(self.dat[:,-1])
            if numpy.any(numpy.logical_and(-isnanpt, -isnandat)):
                self.dat = numpy.append(self.dat, pt, axis = 1)
            else:
                self.dat[isnandat,-1] = pt[isnandat,-1]
        def save(self, fileName):
            numpy.savetxt(fileName, numpy.nan_to_num(self.dat).transpose(), fmt='%.4E')
        def reset(self):
            self.dat = numpy.array([[numpy.NaN],[numpy.NaN],[numpy.NaN]])    #Control, Measured, Sigma
        @property
        def C(self):
            return numpy.nan_to_num(self.dat[0])
        @property
        def M(self):
            return numpy.nan_to_num(self.dat[1])
        @property
        def S(self):
            return numpy.nan_to_num(self.dat[2])

    class Plotter(QtCore.QObject):
        plotSignal = QtCore.pyqtSignal(int)
        def __init__(self, HData):
            QtCore.QObject.__init__(self)
            self.Data = HData
            self.plotSignal.connect(self.updatePlotSlot)
        def setFig(self, figN = 1, xlim = None, ylim = None, xlabel = '', ylabel = '', cpl = 'bo-', errBars = False):
            f = plt.figure(figN)
            self.line, = plt.plot([],[],cpl)
            plt.xlim(xlim)
            plt.ylim(ylim)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.grid(True)
            self.xlim = xlim
            self.ylim = ylim
            self.xlabel = xlabel
            self.ylabel = ylabel
            self.cpl = cpl
            f.canvas.draw()
        def updatePlotSlot(self, figN = 1):
            f = plt.figure(figN)
            if f.axes == []:
                self.setFig(figN, self.xlim, self.ylim, self.xlabel, self.ylabel, self.cpl)
            self.line.set_data(self.Data.C, self.Data.M)
            f.canvas.draw()
        def plot(self, figN = 1):
            self.plotSignal.emit(figN)

    def __init__(self):
        self._lg = self.Logger()
        self.log = self._lg.log
        self.Data = self.HData()
        self._pl = self.Plotter(self.Data)
        self.plot = self._pl.plot
        self.setFig = self._pl.setFig
        self.FC = FieldControler(parent = self)
        self.VC = VSMControler(parent = self)
        
        self.MainThread = self.MThread()
        self.MainThread.start()
        
    def FreqCurve(self, Field = 2000, fcurve = 'auto', file = None):
        aS = self.MainThread.addFunction
        F = lambda f: self.Data.addPoint(C = f)
        A = lambda a: self.Data.addPoint(M = a)
        logTime = lambda t: self.log(t, 'b')
        
        self.Data.reset()
        if fcurve == 'auto':
            fcurve = numpy.linspace(20, 60, 41)
        pts = len(fcurve)

        TC = self.VC.LockIn.TC        
        self.log('\n Frequency Curve: \n', 'r')
        aS(self.FC.setField, 0, Field)
        
        yl = self.VC.LockIn.SEN * self.VC.emu_per_V
        xm = numpy.min(fcurve)
        xM = numpy.max(fcurve)
        figN = 11
        self.setFig(figN, [xm, xM], [0, yl], 'Frequency (Hz)', 'Amplitude')

        for i in range(pts):
            aS(self.log, 0, str(i+1) + '/' + str(pts) + ' --- ', 'b', EOL = '')
            aS(time.ctime, logTime)
            aS(self.VC.confDriver, 0, OscFrec = fcurve[i])
            aS(self.VC.getAmplitude, A, n = 50, iniDelay = 5*TC, measDelay = TC/5.0)
            aS(F, 0, fcurve[i])
            aS(self.plot, 0, figN)
        aS(self.FC.TurnOff)
        
        if file != None:
            aS(self.log, 0, 'Saving data to : ' + file, 'g')
            aS(self.Data.save, 0, file)
        aS(self.log, 0, 'Frequency Curve FINISHED \n', 'r')
        aS(self.FC.BEEP)
            
    def HistCurve(self, MaxField, np = 201, curve = 'auto', file = None, Res = 'Low'):
        aS = self.MainThread.addFunction
        H = lambda h: self.Data.addPoint(C = h)
        MS = lambda ms: self.Data.addPoint(M = ms[0], S = ms[1])
        logTime = lambda t: self.log(t, 'b')
        
        self.Data.reset()
        if curve == 'auto':
            curve = numpy.concatenate((
                            numpy.linspace(MaxField, -MaxField, np/2),
                            numpy.linspace(-MaxField, MaxField, np/2)))
        pts = len(curve)
        
        TC = self.VC.LockIn.TC
        self.log('\n Hysteresis Curve: \n', 'r')
        aS(self.FC.setField, 0, -MaxField)
        aS(self.FC.setField, 0, MaxField)

        yl = self.VC.LockIn.SEN * self.VC.emu_per_V
        xm = numpy.min(curve)
        xM = numpy.max(curve)
        figN = 10
        self.setFig(figN, [xm, xM], [-yl, yl], 'Field (Oe)', 'Magnetization x G (e.m.u.)')
        
        for i in range(pts):
            aS(self.log, 0, str(i+1) + '/' + str(pts) + ' --- ', 'b', EOL = '')
            aS(time.ctime, logTime)
            aS(self.FC.setField, 0, curve[i])

            if Res == 'High':
                aS(self.VC.getMagnetization, MS, n = 30, iniDelay = 5*TC, measDelay = TC/3.0, stat = True, tol = 0.02)
            elif Res == 'Med':
                aS(self.VC.getMagnetization, MS, n = 100, measDelay = TC/10.0)
            else:
                aS(self.VC.getMagnetization, MS, n = 50)
            aS(self.FC.getField, H, Res = 'High')
            aS(self.plot, 0, figN)
        aS(self.FC.TurnOff)
        
        if file != None:
            aS(self.log, 0, 'Saving data to : ' + file, 'g')
            aS(self.Data.save, 0, file)
        aS(self.log, 0, 'Hysteresis Curve FINISHED \n', 'r')
        aS(self.FC.BEEP)
    def STOP(self):
        self.MainThread.clear()
        self.log('\n STOPING!! ', 'r')
        self.MainThread.addFunction(self.FC.TurnOff)
        self.MainThread.addFunction(self.FC.BEEP)
        self.MainThread.addFunction(self.log, 0, 'DONE')
        
def defCurve(pis, pfs, ns):
    x = numpy.array([])
    for i in range(len(ns)):
        x = numpy.concatenate((x, numpy.linspace(pis[i],pfs[i],ns[i])))
    return x
        
#Funciones para log y data
def _empty_log(*args, **kwargs):
    pass
def _empty_data(*args, **kwargs):
    pass


vsm = VSM()
def setField(fld):
    vsm.MainThread.addFunction(vsm.FC.setField, 0, fld)
def confInput(Sen = 1, TC = 0.1, AcGain = '0'):
    vsm.VC.confInput(Sen, TC, AcGain)
def HistCurve(MaxField, np = 201, curve = 'auto', file = None,  Res = 'High'):
    vsm.HistCurve(MaxField, np, curve, file,  Res)
def correction(H):
    ft = numpy.array([8.35734782E-10, 2.72820234E-8]) * vsm.VC.emu_per_V
    return numpy.polyval(ft, H)
def Mc(H, M, HSat):
    Ha1 = H[H >= HSat]
    Ma1 = M[H >= HSat]
    ft1 = numpy.polyfit(Ha1, Ma1, 1)
    Ha2 = H[H <= -HSat]
    Ma2 = M[H <= -HSat]
    ft2 = numpy.polyfit(Ha2, Ma2, 1)
    ft = 0.5*(ft1[0]+ft2[0])
    return M - ft*H
def PlotFile(file, opts = 'bo-', leg = 'Legenda'):
    f = plt.figure(13)
    (h,m) = loadtxt(file, unpack = True)
    f.subplots_adjust(left = 0.2)
    if leg != 'Legenda':
        plt.plot(h, m, opts, label = leg)
    else:
        plt.plot(h, m, opts)
    plt.grid(True)
    plt.xlim(h.min(), h.max())
    plt.xlabel('Field (Oe)')
    plt.ylabel('Magnetization (e.m.u.)')
    if leg != 'Legenda':
        plt.legend(loc = 'lower right')

def HELP():
    print 'Commands'
    print 'setField(field)'
    print 'confInput(Sen = 1E-3, TC = 0.1, AcGain = \'0\')'
    print 'HistCurve(MaxField, np = 201, curve = \'auto\', file = None,  Res = \'Low\')'
    print 'HistCurve(MaxField, curve = crv, file = \'file.txt\',  Res = \'Med\')'
    print 'HistCurve(2000, 0, crv, \'file.txt\', \'High\')'
    print 'defCurve(pis, pfs, ns)'
    print 'PlotFile(file, opts = \'bo-\', leg = \'Legenda\')'
    print '    '
    print 'correction(H)'
    print 'Mc(H, M, HSat)'