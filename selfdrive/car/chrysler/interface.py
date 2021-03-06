#!/usr/bin/env python3
from cereal import car
from selfdrive.controls.lib.drive_helpers import EventTypes as ET, create_event
from selfdrive.car.chrysler.values import Ecu, ECU_FINGERPRINT, CAR, FINGERPRINTS
from selfdrive.car import STD_CARGO_KG, scale_rot_inertia, scale_tire_stiffness, is_ecu_disconnected, gen_empty_fingerprint
from selfdrive.car.interfaces import CarInterfaceBase
from common.params import Params

class CarInterface(CarInterfaceBase):
  @staticmethod
  def compute_gb(accel, speed):
    return float(accel) / 3.0

  @staticmethod
  def get_params(candidate, fingerprint=gen_empty_fingerprint(), has_relay=False, car_fw=[]):
    ret = CarInterfaceBase.get_std_params(candidate, fingerprint, has_relay)
    ret.carName = "chrysler"
    ret.safetyModel = car.CarParams.SafetyModel.chrysler

    # Chrysler port is a community feature, since we don't own one to test
    ret.communityFeature = True

    # Speed conversion:              20, 45 mph
    ret.wheelbase = 3.089  # in meters for Pacifica Hybrid 2017
    ret.steerRatio = 16.2 # Pacifica Hybrid 2017
    ret.mass = 1964. + STD_CARGO_KG  # kg curb weight Pacifica General
    ret.lateralTuning.pid.kpBP, ret.lateralTuning.pid.kiBP = [[9., 20.], [9., 20.]]
    ret.lateralTuning.pid.kpV, ret.lateralTuning.pid.kiV = [[0.15,0.30], [0.03,0.05]]
    ret.lateralTuning.pid.kf = 0.00006   # full torque for 10 deg at 80mph means 0.00007818594
    ret.steerActuatorDelay = 0.1
    ret.steerRateCost = 0.7
    ret.steerLimitTimer = 0.4

    if candidate in (CAR.JEEP_CHEROKEE, CAR.JEEP_CHEROKEE_2019):
      ret.wheelbase = 2.91  # in meters
      ret.steerRatio = 12.7
      ret.steerActuatorDelay = 0.2  # in seconds

    if candidate in (CAR.CHRYSLER_300_2018):
      ret.wheelbase = 3.05308 # in meters
      ret.steerRatio = 15.5 # 2013 V-6 (RWD) — 15.5:1 V-6 (AWD) — 16.5:1 V-8 (RWD) — 15.5:1 V-8 (AWD) — 16.5:1
      ret.mass = 1828.0 + STD_CARGO_KG # 2013 V-6 RWD
      ret.lateralTuning.pid.kf = 0.00006   # full torque for 10 deg at 80mph means 0.00007818594
      ret.steerLimitTimer = 0.1

      
  # "AndrewSteerRateCost": [TxType.PERSISTENT],
  # "AndrewSteerLimitTimer": [TxType.PERSISTENT],
  # "AndrewINDIInnerLoopGain": [TxType.PERSISTENT],
  # "AndrewINDIOuterLoopGain": [TxType.PERSISTENT],
  # "AndrewINDIActuatorEffectiveness": [TxType.PERSISTENT],
  # "AndrewINDITimeConstant": [TxType.PERSISTENT],

    # # Change to try on all vehicles!!
    # params = Params()

    # if params.get("AndrewSteerRateCost") is None:
    #   params.put("AndrewSteerRateCost", "1.0")
    
    # if params.get("AndrewSteerLimitTimer") is None:
    #   params.put("AndrewSteerLimitTimer", "0.8")
    
    # if params.get("AndrewINDIInnerLoopGain") is None:
    #   params.put("AndrewINDIInnerLoopGain", "2.53")
    
    # if params.get("AndrewINDIOuterLoopGain") is None:
    #   params.put("AndrewINDIOuterLoopGain", "0.92")
    
    # if params.get("AndrewINDIActuatorEffectiveness") is None:
    #   params.put("AndrewINDIActuatorEffectiveness", "1.35")
    
    # if params.get("AndrewINDITimeConstant") is None:
    #   params.put("AndrewINDITimeConstant", "1.0")

    # if params.get("AndrewSteerActuatorDelay") is None:
    #   params.put("AndrewSteerActuatorDelay", "0.1")
    
    # ret.steerActuatorDelay =  float(params.get("AndrewSteerActuatorDelay", encoding='utf8'))

    # ret.steerRateCost = float(params.get("AndrewSteerRateCost", encoding='utf8'))
    # ret.steerLimitTimer = float(params.get("AndrewSteerLimitTimer", encoding='utf8'))

    # ret.lateralTuning.init('indi')
    # ret.lateralTuning.indi.innerLoopGain = float(params.get("AndrewINDIInnerLoopGain", encoding='utf8'))
    # ret.lateralTuning.indi.outerLoopGain = float(params.get("AndrewINDIOuterLoopGain", encoding='utf8'))
    # ret.lateralTuning.indi.timeConstant = float(params.get("AndrewINDITimeConstant", encoding='utf8'))
    # ret.lateralTuning.indi.actuatorEffectiveness = float(params.get("AndrewINDIActuatorEffectiveness", encoding='utf8'))
    

    ret.centerToFront = ret.wheelbase * 0.44

    ret.minSteerSpeed = 3.8  # m/s
    if candidate in (CAR.PACIFICA_2019_HYBRID, CAR.PACIFICA_2020, CAR.JEEP_CHEROKEE_2019):
      # TODO allow 2019 cars to steer down to 13 m/s if already engaged.
      ret.minSteerSpeed = 17.5  # m/s 17 on the way up, 13 on the way down once engaged.

    # starting with reasonable value for civic and scaling by mass and wheelbase
    ret.rotationalInertia = scale_rot_inertia(ret.mass, ret.wheelbase)

    # TODO: start from empirically derived lateral slip stiffness for the civic and scale by
    # mass and CG position, so all cars will have approximately similar dyn behaviors
    ret.tireStiffnessFront, ret.tireStiffnessRear = scale_tire_stiffness(ret.mass, ret.wheelbase, ret.centerToFront)

    ret.enableCamera = is_ecu_disconnected(fingerprint[0], FINGERPRINTS, ECU_FINGERPRINT, candidate, Ecu.fwdCamera) or has_relay

    return ret

  # returns a car.CarState
  def update(self, c, can_strings):
    self.dp_load_params('chrysler')
    # ******************* do can recv *******************
    self.cp.update_strings(can_strings)
    self.cp_cam.update_strings(can_strings)

    ret = self.CS.update(self.cp, self.cp_cam)

    ret.canValid = self.cp.can_valid and self.cp_cam.can_valid

    # speeds
    ret.steeringRateLimited = self.CC.steer_rate_limited if self.CC is not None else False

    ret.buttonEvents = []

    # events
    events = self.create_common_events(ret, extra_gears=[car.CarState.GearShifter.low], gas_resume_speed=2.)

    if ret.vEgo < self.CP.minSteerSpeed:
      events.append(create_event('belowSteerSpeed', [ET.WARNING]))

    ret.events = events

    # copy back carState packet to CS
    self.CS.out = ret.as_reader()

    return self.CS.out

  # pass in a car.CarControl
  # to be called @ 100hz
  def apply(self, c):

    if (self.CS.frame == -1):
      return [] # if we haven't seen a frame 220, then do not update.

    can_sends = self.CC.update(c.enabled, self.CS, c.actuators, c.cruiseControl.cancel, c.hudControl.visualAlert)

    return can_sends
