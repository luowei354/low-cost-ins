%% LCNAVI Inertial Navigation class %
%
% LCNAVI class is intended for low-cost aided-inertial navigation
% applications. In that sense, we aim for low precision sensors commonly
% found in micro-air vehicles (nevertheless, we still pursue great accuracy
% and some statistical consistency). Therefore, the navigation routines
% found herein does not account for Earth rotation and assume small
% displacements. The complete science and motivation behind this algorithm 
% can be found in
%
% [!] L.R. Lustosa, S. Pizziol, F. Defay, J.M. Moschetta, "An error model
% of a complementary filter for use in Bayesian filtering - the CF-EKF
% filter", 2016
% (to be published)
%
% Currently supported sensors: accelerometers, rate-gyros, magnetometers, 
% GNSS receivers (loosely coupled).
%
classdef lcNavi < handle
    
    %% Navigator internal constants
    properties (Constant, GetAccess='public')
        % default magnetics based on [www.magnetic-declination.com]
        MAG_INC_ISAE = 0*(59.0+5.0/60.0)*pi/180.0; % inclination at ISAE (rad)
        MAG_DEC_ISAE = 0*(18.0/60.0)*pi/180.0; % declination at ISAE (rad)
        MAG_ISAE = 4.62329e-5; % magnetics magnitude at ISAE (T)
        % default location based on [www.magnetic-declination.com]
        LLA_ISAE = [ (43+33/60+57.4/60^2)*pi/180; ... % [rad] lat
                     ( 1+28/60+29.3/60^2)*pi/180; ... % [rad] lon
                     0 ]; % [m] altitude
        % default gravity value based on [http://www.physicsclassroom.com/class/1DKin/Lesson-5/Acceleration-of-Gravity]
        G = 9.80653; % down component at Toulouse, France.
        % default IMU sampling rate
        %% Complementary filter parameters default
        CF_GAIN_ACC = 0.4/9.80653;
        CF_GAIN_MAG = 0.8/4.62329e-5;
        %% orientation definition flags
        MODE_QUATERNION = 1; % 4x1 quaternion description
        MODE_DCM = 2; % 3x3 dcm description
        MODE_EULER = 3; % 3x1 [roll; pitch; yaw] in 3-2-1 order euler
        %% INS operation modes
        % in OP_IDLE mode, the IMU disregards any data it receives and do
        % nothing. This mode can be used to process initial meaningless
        % data that are often found in beginning of flight tests.
        OP_IDLE = 0;
        % in OP_DYNAMIC_CALIBRATION mode, the IMU is expected to realize as many different
        % attitude configurations as possible as to the magnetometer to
        % infer mins and maxs of local magnetic field. This additionally
        % cancels out the drone magnetic residues on the magnetometer.
        OP_DYNAMIC_CALIBRATION = 1;
        % in OP_STATIC_CALIBRATION (formerly OP_GYR_CALIBRATION) mode, the IMU is expected to be still with
        % respect to Earth and perform rate-gyro bias estimation. All other
        % sensor noise covariances can be estimated at this moment as well
        % (for Kalman filtering purposes).
        OP_STATIC_CALIBRATION = 2;
        % in OP_CRUISE, the IMU has no a priori information over its
        % movement and compute navigation variables using a EKF
        % over a complementary filter by means of sensor measurements.
        OP_CRUISE = 3;
        % in OP_ALIGNMENT, the... *** TODO *** documentation
        OP_ALIGNMENT = 4;
        %% Kalman Filter operation modes
        % the following mode defines the naive implementation of the Kalman
        % filter equations, which could be not suited to embedded
        % applications. This one uses the Pk matrix.
        % THIS IS THE ONLY ONE IMPLEMENTED TODAY!
        KALMAN_NAIVE = 1;
        % this mode is the more robust implementation (UD-factorization).
        % This mode will use the Uk and Dk matrices.
        % THIS IS NOT YET IMPLEMENTED! COME BACK LATER!
        KALMAN_UD = 2;
        %% suboptimal vibration levels expected in vehicle trajectory
        % this value can be obtained by evaluating a sample trajectory and
        % gradually increasing VIBRATION_SIG such that predictions (no
        % additional sensor updates!) yield a consistent filter.
        VIBRATION_SIG = 2.0;
    end
    
    %% Navigator internal variables
    properties (GetAccess='private', SetAccess='private')
        %% INS State variables
        % Kinematics
        Pos = [0;0;0]; % position (pn,pe,pd) in local NED (m)
        Vel = [0;0;0]; % velocity wrt Earth in local NED (m/s)
        Quat = [1;0;0;0]; % quaternion of body wrt local NED
        % INS operation mode (more info on constants definitions)
        OpMode = lcNavi.OP_IDLE; 
        %% Kalman filter parameters
        % each kalman filter reset yields the following covariance (value does not change in time)
        % ps: this variable is accessed directly in the code (i.e., no set and get functions)
        P0 = diag([4/3*2.5*ones(3,1); 0.3*ones(3,1); 5*pi/180*ones(3,1); 0.03*ones(3,1); (0.2*pi/180)*ones(3,1); 0.005*ones(3,1) ])^2;
        % the kalman filter covariance at any given moment is the follows
        % (values changes in time)
        Pk = diag([0.1*ones(3,1); 0.1*ones(3,1); 0.01*ones(3,1); 0.05*ones(3,1); 1e-4*ones(3,1); 1e-6*ones(3,1) ])^2;
        % we define x = [dp; dv; psi; bias_acc; bias_gyr; bias_mag] (loosely coupled, simplified)
        % we define w = [noise_acc; noise_gyr; noise_mag; noise_vibrations]
        % for more info on the following formulas, check my article
        xk = zeros(18,1); % We have no idea over the estimates, so we set zero
        KalmanMode = lcNavi.KALMAN_NAIVE; % inits in naive mode
        % sensor calibration should be done or not? (i.e., KalmanCalibration = 1 or 0)
        % ps: this variable is accessed directly in the code (i.e., no set and get functions)
        KalmanCalibration = 0; 
        % here comes a lot of default sensor statistics
        Qk = diag([zeros(1,3) zeros(1,3) zeros(1,3) zeros(1,3) ]);
        Rk_pos = (4/3*2.5)^2*eye(3);
        Rk_vel = (4/3*0.1)^2*eye(3);
        %% Bias parameters
        % these values compensate sensors bias directly. notice that this
        % is not the EKF estimates! however, EKF calibration is performed
        % by transporting the values of the EKF to the following variables.
        % INS stand-alone algorithm takes into account these values,
        % whether or not an EKF is present.
        gyr_Bias = zeros(3,1); % [rad/s]
        mag_Bias = zeros(3,1); % [T]
        acc_Bias = zeros(3,1); % [m/s^2]
        N_gyr_meas = 0;
        N_mag_meas = 0;
        N_acc_meas = 0;
        maxMag = zeros(3,1); % for dynamic calibration x y z in body
        minMag = zeros(3,1); % for dynamic calibration x y z in body
        %% Position of local NED definition (at ISAE)
        LLA0 = lcNavi.LLA_ISAE; % LLA of local NED frame definition [rad,rad,m]
        %% Local magnetic field parameters
        % initialization assumes drone is in ISAE, TOULOUSE.
        % this should be reconfigured depending on the location of the
        % drone.
        inclination = lcNavi.MAG_INC_ISAE; % [rad]
        declination = lcNavi.MAG_DEC_ISAE; % [rad] 
        mag_magnitu = lcNavi.MAG_ISAE;     % [ T ] 
        %% IMU sampling rate and time
        tempo = -1; % current timestamp in [sec]
        %% Complementary filter parameters
        CF_Kg = lcNavi.CF_GAIN_ACC; % default values init
        CF_Km = lcNavi.CF_GAIN_MAG; % default values init
    end
    
    %% Navigator public methods
    methods (Access=public)
        
        function predict(obj, gyr, acc, mag, timestamp)
            %% PREDICT: update INS based on IMU and magnetometer (normal operation)
            %
            % inputs: gyr - 3x1 Real - uncompensated gyro measurement [rad/s]
            %         acc - 3x1 Real - uncompensated accelerometer measurement [m/s^2]
            %         mag - 3x1 Real - uncompensated magnetometer measurement [T]
            %         timestamp - Scalar Real - timestamp [sec]
            
            % let's make sure we are in right operation mode
            obj.setOpMode(obj.OP_CRUISE);
            % compute difference in time since last call to INS functions
            dt = timestamp - obj.getTempo();
            % for starters, update complementary filter
            obj.updateCF(gyr, acc, mag, dt);
            % compute error model for Kalman filter usage
            [Fk, Gk] = obj.getINSErrorModel(gyr, acc, mag, dt);
            % get current statistics
            x = obj.getKalmanState();
            P = obj.getKalmanCov();
            % get process noise statistics
            Q = obj.getKalmanProcNoise()*dt^2; 
            % update statistics using Kalman prediction formulas
            x_p = Fk*x;
            P_p = Fk*P*Fk' + Gk*Q*Gk';
            obj.setKalmanState(x_p);
            obj.setKalmanCov(P_p);
            % update timestamp
            obj.setTempo(timestamp);
        end
        
        function looselyCoupledPosFix(obj, lat, lon, h)
            %% LOOSELYCOUPLEDPOSFIX: GNSS position fix
            %
            % inputs: lat - Real scalar - GNSS WGS-84 latitude [rad]
            %         lon - Real scalar - GNSS WGS-84 longitude [rad]
            %         alt - Real scalar - GNSS WGS-84 altitude [m]
            
            %% let's make sure we are in right operation mode
            obj.setOpMode(obj.OP_CRUISE);
            %% first conversion of WGS-84 to local NED
            Pos_GNSS = obj.lla2posNED([lat lon h]');
            %% then Kalman filter update formula
            Hk = obj.getGNSSlooselyPosModel();
            P = obj.getKalmanCov();
            x = obj.getKalmanState();
            % construct sensor output
            yk = obj.getPosition() - Pos_GNSS;
            % get sensor covariance matrix
            Rk = obj.getKalmanPosNoise();
            % Kalman formulas
            res = yk - Hk*x;
            S = Hk*P*Hk' + Rk;
            x_updated = x + P*Hk'*(S\res); 
            P_updated = (eye(18) - P*Hk'*(S\Hk))*P;
            %% update estimated state and covariance
            obj.setPosition( obj.getPosition() - x_updated(1:3) );
            x_updated(1:3) = zeros(3,1);
            obj.setVelocity( obj.getVelocity() - x_updated(4:6) );
            x_updated(4:6) = zeros(3,1);
            obj.setAttitude( obj.MODE_DCM, angle2dcm(x_updated(9),x_updated(8),x_updated(7))*obj.getAttitude(obj.MODE_DCM)  );
            x_updated(7:9) = zeros(3,1);
            if obj.KalmanCalibration
                obj.setAccEstBias( obj.getAccEstBias() + x_updated(10:12) );
                x_updated(10:12) = zeros(3,1);
                obj.setGyroEstBias( obj.getGyroEstBias() + x_updated(13:15) );
                x_updated(13:15) = zeros(3,1);
                obj.setMagEstBias( obj.getMagEstBias() + x_updated(16:18) );
                x_updated(16:18) = zeros(3,1);
            end
            obj.setKalmanState(x_updated);
            obj.setKalmanCov(P_updated);
        end
        
        function looselyCoupledVelFix(obj, vn, ve, vd)
            %% LOOSELYCOUPLEDVELFIX: GNSS velocity fix
            %
            % inputs: vn - Real scalar - GNSS velocity north [m/s]
            %         ve - Real scalar - GNSS velocity east  [m/s]
            %         vd - Real scalar - GNSS velocity down  [m/s]
        
            %% let's make sure we are in right operation mode
            setOpMode(obj, obj.OP_CRUISE);
            %% Kalman filter update formula
            Hk = getGNSSlooselyVelModel(obj);
            P = obj.getKalmanCov();
            x = obj.getKalmanState();
            Vel_GNSS = [vn; ve; vd];
            % construct sensor output
            yk = obj.getVelocity() - Vel_GNSS;
            % get sensor covariance matrix
            Rk = obj.getKalmanVelNoise();
            % Kalman formulas
            res = yk - Hk*x;
            S = Hk*P*Hk' + Rk;
            x_updated = x + P*Hk'*(S\res);
            P_updated = (eye(18,18) - P*Hk'*(S\Hk))*P;
            %% performs CF correction and CF sensors calibration
            % correction
            obj.setPosition( obj.getPosition() - x_updated(1:3) );
            x_updated(1:3) = zeros(3,1);
            obj.setVelocity( obj.getVelocity() - x_updated(4:6) );
            x_updated(4:6) = zeros(3,1);
            obj.setAttitude( obj.MODE_DCM, angle2dcm(x_updated(9),x_updated(8),x_updated(7))*obj.getAttitude(obj.MODE_DCM)  );
            x_updated(7:9) = zeros(3,1);
            % calibration
            if obj.KalmanCalibration
                obj.setAccEstBias( obj.getAccEstBias() + x_updated(10:12) );
                x_updated(10:12) = zeros(3,1);
                obj.setGyroEstBias( obj.getGyroEstBias() + x_updated(13:15) );
                x_updated(13:15) = zeros(3,1);
                obj.setMagEstBias( obj.getMagEstBias() + x_updated(16:18) );
                x_updated(16:18) = zeros(3,1);
            end
            %% update estimated state and covariance
            obj.setKalmanState(x_updated);
            obj.setKalmanCov(P_updated);
        end
        
        function resetEKF(obj)
            %% RESETEKF: reset Kalman filter covariance matrix
            %
            % this should be used in case statistical inconsistency is
            % detected by the user. another approach is to periodically
            % reset the filter to make sure it does not get too much
            % optimistic (remember we implement a suboptimal solution here,
            % if the covariance gets too small, we might incur into
            % statistical inconsistency).
            
            obj.Pk = obj.P0;
        end
        
        function staticMagCalibration(obj, mag)
            %% STATICMAGCALIBRATION: calibrate sensors when INS is not moving
            %
            % NOTICE that one must perform the dynamic calibration again after
            % a static calibration! The static calibration ruins the dynamic
            % calibration!
            %
            % inputs: mag - 3x1 Real - uncompensated magnetometer measurement [T]
        
            % is it the first point in the calibration procedure?
            if obj.getOpMode ~= obj.OP_STATIC_CALIBRATION
                % reset calibration values
                bias=zeros(3,1);
                obj.setGyroEstBias(bias);
                obj.setAccEstBias(bias);
                obj.setMagEstBias(bias);
                obj.N_gyr_meas = 0;
                obj.N_acc_meas = 0;
                obj.N_mag_meas = 0;
                obj.Qk = diag([zeros(9,1); obj.VIBRATION_SIG^2*ones(3,1)]);
                % change mode
                obj.setOpMode(obj.OP_STATIC_CALIBRATION);
            end
            % now we do calibration
            N = obj.N_mag_meas;
            prev_mean = obj.getMagEstBias();
            prev_mean_sum = N*prev_mean;
            prev_cov = obj.Qk(7:9,7:9); % magnetometer position
            prev_cov_sum = N*prev_cov;
            N = N + 1;
            new_mean_sum = prev_mean_sum + mag;
            new_mean = new_mean_sum/N;
            new_cov_sum = prev_cov_sum + (mag-new_mean)*(mag-new_mean)';
            new_cov = new_cov_sum/N;
            % update estimates
            obj.N_mag_meas = N;
            obj.setMagEstBias(new_mean);
            obj.Qk(7:9,7:9) = new_cov;
        end
        
        function staticGyroCalibration(obj, gyr)
            %% STATICGYROCALIBRATION: calibrate sensors when INS is not moving
            %
            % NOTICE that one must perform the dynamic calibration again after
            % a static calibration! The static calibration ruins the dynamic
            % calibration! (except for rategyro bias, it will be ok already)
            %
            % inputs: gyr - 3x1 Real - uncompensated gyro measurement [rad/s]
            
            % is it the first point in the calibration procedure?
            if obj.getOpMode ~= obj.OP_STATIC_CALIBRATION
                % reset calibration values
                obj.setGyroEstBias(zeros(3,1));
                obj.setAccEstBias(zeros(3,1));
                obj.setMagEstBias(zeros(3,1));
                obj.N_gyr_meas = 0;
                obj.N_acc_meas = 0;
                obj.N_mag_meas = 0;
                obj.Qk = diag([zeros(9,1); obj.VIBRATION_SIG^2*ones(3,1)]);
                % change mode
                obj.setOpMode(obj.OP_STATIC_CALIBRATION);
            end
            % now we do calibration
            N = obj.N_gyr_meas;
            prev_mean = obj.getGyroEstBias();
            prev_mean_sum = N*prev_mean;
            prev_cov = obj.Qk(4:6,4:6); % gyro position
            prev_cov_sum = N*prev_cov;
            N = N + 1;
            new_mean_sum = prev_mean_sum + gyr;
            new_mean = new_mean_sum/N;
            new_cov_sum = prev_cov_sum + (gyr-new_mean)*(gyr-new_mean)';
            new_cov = new_cov_sum/N;
            % update estimates
            obj.N_gyr_meas = N;
            obj.setGyroEstBias(new_mean);
            obj.Qk(4:6,4:6) = new_cov;
        end
        
        function staticAccCalibration(obj, acc)
            %% STATICACCCALIBRATION: calibrate sensors when INS is not moving
            %
            % NOTICE that one must perform the dynamic calibration again after
            % a static calibration! The static calibration ruins the dynamic
            % calibration! (except for rategyro bias, it will be ok already)
            %
            % inputs: acc - 3x1 Real - uncompensated accelerometer measurement [m/s^2]
            
            % is it the first point in the calibration procedure?
            if obj.getOpMode ~= obj.OP_STATIC_CALIBRATION
                % reset calibration values
                obj.setGyroEstBias(zeros(3,1));
                obj.setAccEstBias(zeros(3,1));
                obj.setMagEstBias(zeros(3,1));
                obj.N_gyr_meas = 0;
                obj.N_acc_meas = 0;
                obj.N_mag_meas = 0;
                obj.Qk = diag([zeros(9,1); obj.VIBRATION_SIG^2*ones(3,1)]);
                % change mode
                obj.setOpMode(obj.OP_STATIC_CALIBRATION);
            end
            % now we do calibration
            N = obj.N_acc_meas;
            prev_mean = obj.getAccEstBias();
            prev_mean_sum = N*prev_mean;
            prev_cov = obj.Qk(1:3,1:3); % acc position
            prev_cov_sum = N*prev_cov;
            N = N + 1;
            new_mean_sum = prev_mean_sum + acc;
            new_mean = new_mean_sum/N;
            new_cov_sum = prev_cov_sum + (acc-new_mean)*(acc-new_mean)';
            new_cov = new_cov_sum/N;
            % update estimates
            obj.N_acc_meas = N;
            obj.setAccEstBias(new_mean);
            obj.Qk(1:3,1:3) = new_cov;
        end
        
        function dynamicMagCalibration(obj, mag)
            %% DYNAMICMAGCALIBRATION: calibrate sensors when INS is moving
            %
            % NOTICE that one must perform the dynamic calibration again after
            % a static calibration! The static calibration ruins the dynamic
            % calibration! (except for rategyro bias, it will be ok already)
            %
            % inputs: mag - 3x1 Real - uncompensated magnetometer measurement [T]
            
            % is it the first point in the dynamic calibration procedure?
            if obj.getOpMode ~= obj.OP_DYNAMIC_CALIBRATION
                % reset dynamic calibration values
                obj.setAccEstBias(zeros(3,1)); % this fix the static calib
                obj.setMagEstBias(zeros(3,1)); % this fix the static calib
                obj.maxMag = -1e10*ones(3,1);
                obj.minMag = 1e10*ones(3,1);
                % change mode
                obj.setOpMode(obj.OP_DYNAMIC_CALIBRATION);
            end
            
            % check all XYZ magnetometer channels
            for i = 1:3
                if mag(i) > obj.maxMag(i)
                    obj.maxMag(i) = mag(i);
                end
                if mag(i) < obj.minMag(i)
                    obj.minMag(i) = mag(i);
                end
            end
            
            % update magnetometer bias estimate
            bias = (obj.minMag+obj.maxMag)/2;
            obj.setMagEstBias(bias);
            
        end
        
        function initialAlignment(obj, gyr, acc, mag, timestamp)
            %% INITIALALIGNMENT: search for INS initial attitude
            %
            % NOTICE: this must be done always, and after static and dynamic
            % calibration procedures, and before cruise flight, until the
            % attitude reach convergence! for better results, perform
            % calibration with the drone stationary.
            %
            % inputs: gyr - 3x1 Real - uncompensated gyro measurement [rad/s]
            %         acc - 3x1 Real - uncompensated accelerometer measurement [m/s^2]
            %         mag - 3x1 Real - uncompensated magnetometer measurement [T]
            %         timestamp - Scalar Real - time [sec]
            
            
            if obj.getOpMode() ~= obj.OP_ALIGNMENT
                % change mode 
                obj.setOpMode(obj.OP_ALIGNMENT);
                % from here on, time flows normally and need to be iniatilly
                % set for correct intial navigation increment of time
                obj.setTempo(timestamp);
            end
            % compute difference in time since last call to INS functions
            dt = timestamp - obj.getTempo();
            % for starters, update complementary filter!
            % this will yield crazy position and velocity measurements that
            % need to be corrected setting the local NED LLA by means of
            % GNSS in the user code!
            obj.updateCF(gyr, acc, mag, dt);
            obj.setPosition(zeros(3,1)); % just alignement here!!!
            obj.setVelocity(zeros(3,1)); % just alignement here!!!
            % updates current time
            obj.setTempo(timestamp);
            
        end
        
        function xk = getKalmanState(obj)
            %% GETKALMANSTATE: gets Kalman estimation of state
            %
            % outputs: xk - 18x1 Real - current estimate of state [SI units]
            
            % we define x = [dp; dv; psi; bias_acc; bias_gyr; bias_mag] (loosely coupled, simplified)
            % we define w = [noise_acc; noise_gyr; noise_mag; noise_vibrations]
            xk = obj.xk;
        end
        
        function setKalmanState(obj, xk)
            %% SETKALMANSTATE: sets Kalman estimation of state
            %
            % inputs: xk - 18x1 Real - current estimate of state [SI units]
            
            % we define x = [dp; dv; psi; bias_acc; bias_gyr; bias_mag] (loosely coupled, simplified)
            % we define w = [noise_acc; noise_gyr; noise_mag; noise_vibrations]
            obj.xk = xk;
        end
        
        function Pk = getKalmanCov(obj)
            %% GETKALMANCOV: gets Kalman estimation covariance
            %
            % outputs: Pk - 18x18 Real - current estimate covariance [SI units]
            
            % we define x = [dp; dv; psi; bias_acc; bias_gyr; bias_mag] (loosely coupled, simplified)
            % we define w = [noise_acc; noise_gyr; noise_mag; noise_vibrations]
            Pk = obj.Pk;
        end
        
        function setKalmanCov(obj, Pk)
            %% SETKALMANCOV: sets Kalman estimation covariance
            %
            %  inputs: Pk - 18x18 Real - current estimate covariance [SI units]
            
            % we define x = [dp; dv; psi; bias_acc; bias_gyr; bias_mag] (loosely coupled, simplified)
            % we define w = [noise_acc; noise_gyr; noise_mag; noise_vibrations]
            obj.Pk = Pk;
        end
        
        function Rk = getKalmanPosNoise(obj)
            %% GETKALMANPOSNOISE: gets Kalman GNSS position noise matrix
            %
            % outputs: Rk - 3x3 - position noise covariance [SI units]
            Rk = obj.Rk_pos;
        end
        
        function Rk = getKalmanVelNoise(obj)
            %% GETKALMANVELNOISE: gets Kalman GNSS velocity noise matrix
            %
            % outputs: Rk - 3x3 - velocity noise covariance [SI units]
            Rk = obj.Rk_vel;
        end
        
        function Qk = getKalmanProcNoise(obj)
            %% GETKALMANPROCNOISE: gets Kalman process noise matrix
            %
            % outputs: Qk - 12x12 - process noise covariance [SI units]
            
            % we define x = [dp; dv; psi; bias_acc; bias_gyr; bias_mag] (loosely coupled, simplified)
            % we define w = [noise_acc; noise_gyr; noise_mag; noise_vibrations]
            Qk = obj.Qk;
        end
        
        function setKalmanProcNoise(obj, Qk)
            %% SETKALMANPROCNOISE: gets Kalman process noise matrix
            %
            % inputs: Qk - 12x12 - process noise covariance [SI units]
            
            % we define x = [dp; dv; psi; bias_acc; bias_gyr; bias_mag] (loosely coupled, simplified)
            % we define w = [noise_acc; noise_gyr; noise_mag; noise_vibrations]
            obj.Qk = Qk;
        end
        
        function setAttitude(obj, dMode, att)
            %% SETATTITUDE: sets orientation of body wrt NED
            %
            % inputs: dMode - Real scalar - attitude description
            %         att - can be either 4x1, 3x3, 3x1 [rad]
            % NOTICE: dMode \in {MODE_QUATERNION, MODE_DCM, MODE_EULER}
            % In euler angles: 3x1 [roll; pitch; yaw] in 3-2-1 order
        
            % quaternion mode: input att is 4x1 Real
            if dMode == lcNavi.MODE_QUATERNION
                obj.Quat = att;
                return
            end
            % dcm mode mode: input att is 3x3 Real
            if dMode == lcNavi.MODE_DCM
                obj.Quat = dcm2quat(att)';
                return
            end
            % euler angles mode: input att is 3x1 Real
            if dMode == lcNavi.MODE_EULER
                obj.Quat = angle2quat(att(3), att(2), att(1))';
                return
            end
            % the code was not supposed to reach here!!
            error('Invalid attitude mode description.')
        end
        
        function att = getAttitude(obj, dMode)
            %% GETATTITUDE: gets orientation of body wrt NED
            %
            % inputs:  dMode - Real scalar - attitude description
            % outputs: att - can be either 4x1, 3x3, 3x1 [rad]
            % NOTICE: dMode \in {MODE_QUATERNION, MODE_DCM, MODE_EULER}
            % In euler angles: 3x1 [roll; pitch; yaw] in 3-2-1 order
        
            % quaternion mode: input att is 4x1 Real
            if dMode == lcNavi.MODE_QUATERNION
                att = obj.Quat;
                return
            end
            % dcm mode mode: input att is 3x3 Real
            if dMode == lcNavi.MODE_DCM
                att = quat2dcm(obj.Quat');
                return
            end
            % euler angles mode: input att is 3x1 Real
            if dMode == lcNavi.MODE_EULER
                % MATLAB order of angles in vector is different
                att2 = zeros(3,1);
                [att2(1), att2(2), att2(3) ] = quat2angle(obj.Quat');
                % finally translate
                att = [att2(3); att2(2); att2(1)]; 
                return
            end
            % the code was not supposed to reach here!!
            error('Invalid attitude mode description.')
        end
        
        function setVelocity(obj, vel)
            %% SETVELOCITY: sets velocity wrt Earth in NED coords
            %
            % inputs: vel - 3x1 Real - NED velocity [m/s]
            obj.Vel = [vel(1); vel(2); vel(3)];
        end
        
        function vel = getVelocity(obj)
            %% GETVELOCITY: gets velocity wrt Earth in NED coords
            %
            % outputs: vel - 3x1 Real - NED velocity [m/s]
            vel = obj.Vel;
        end
        
        function setPosition(obj, pos)
            %% SETPOSITION: sets position wrt LL0 in NED coords
            %
            % inputs: pos - 3x1 Real - NED position [m]
            obj.Pos = [pos(1); pos(2); pos(3)];
        end

        function pos = getPosition(obj)
            %% GETPOSITION: gets position wrt LL0 in NED coords
            %
            % outputs: pos - 3x1 Real - NED position [m]
            pos = obj.Pos;
        end
        
        function lla = getLLA(obj)
            %% GETLLA: gets drone LLA coordinates
            %
            % this function performs the computation to transform reference
            % LLA0 coordinates and differential NED position into drone LLA
            % position. this function is necessary because we keep position
            % information as NED position vector and eventually LLA
            % coordinates are required instead. NOTICE: the LLAO component
            % in this class is a reference point for the zero-position of
            % the NED position cordinates, and NOT the LLA coordinate of
            % the drone.
            %
            % outputs : lla - 3x1 Real - WGS-84 LLA coordinates [rad,rad,m]
        
            % get actual position
            [lat0, lon0, h0] = obj.getLLA0(); 
            pos_ned = obj.getPosition();
            pn = pos_ned(1);
            pe = pos_ned(2);
            pd = pos_ned(3);
            % transform the NED differentials IN LLA differentials
            a = 6378137; % semi-major axis of Earth (in meters) 
            ec2 = 6.69437999014e-3; % first eccentricity squared
            ec = sqrt(ec2); % first eccentricity
            RN = a*(1-ec^2)/sqrt((1-(ec*sin(lat0))^2)^3); 
            RE = a/sqrt(1-(ec*sin(lat0))^2);
            dlat = pn/(RN + h0);
            dlon = pe/((RE + h0)*cos(lat0));
            dh = -pd;
            % update using differentials in LLA
            lat = dlat + lat0;
            lon = dlon + lon0;
            h   = dh + h0;
            lla = [lat; lon; h];
        end
        
        function setLocalMag(obj, inc, dec, amp)
            %% SETLOCALMAG: sets local magnetic configuration
            %
            % inputs: inc - Real scalar - inclination [rad]
            %         dec - Real scalar - declination [rad]
            %         amp - Real scalar - amplitude [T]
            % NOTICE: the default is ISAE mag field at 2015
            obj.inclination = inc;
            obj.declination = dec;
            obj.mag_magnitu = amp;
        end
        
        function [inc, dec, amp] = getLocalMag(obj)
            %% GETLOCALMAG: gets local magnetic configuration
            %
            % outputs: inc - Real scalar - inclination [rad]
            %          dec - Real scalar - declination [rad]
            %          amp - Real scalar - amplitude [T]
            % NOTICE: the default is ISAE, Toulouse mag field at 2015
            inc = obj.inclination;
            dec = obj.declination;
            amp = obj.mag_magnitu;
        end
        
        function bias = getGyroEstBias(obj) 
            %% GETGYROESTBIAS: gets rate-gyro estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % outputs: bias - 3x1 Real - estimated Bias in Body axes [rad/s]
            bias = obj.gyr_Bias;
        end
        
        function setGyroEstBias(obj, bias) 
            %% SETGYROESTBIAS: sets rate-gyro estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % inputs: bias - 3x1 Real - estimated Bias in Body axes [rad/s]
            obj.gyr_Bias = bias;
        end
        
        function addGyroEstBias(obj, dbias) 
            %% ADDGYROESTBIAS: adds value to rate-gyro estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % inputs: incremental dbias - 3x1 Real - incremental Bias [rad/s]
            obj.gyr_Bias = obj.gyr_Bias + dbias;
        end
        
        function bias = getMagEstBias(obj) 
            %% GETMAGESTBIAS: gets magnetometer estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % outputs: bias - 3x1 Real - estimated Bias in Body axes [T]
            bias = obj.mag_Bias;
        end
        
        function setMagEstBias(obj, bias) 
            %% SETMAGESTBIAS: sets magnetometer estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % inputs: bias - 3x1 Real - estimated Bias in Body axes [T]
            obj.mag_Bias = bias;
        end
        
        function addMagEstBias(obj, dbias) 
            %% ADDMAGESTBIAS: adds value to magnetometer estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % inputs: incremental dbias - 3x1 Real - incremental Bias [T]
            obj.mag_Bias = obj.mag_Bias + dbias;
        end
        
        function bias = getAccEstBias(obj) 
            %% GETACCESTBIAS: gets accelerometer non-EKF estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % outputs: bias - 3x1 Real - estimated Bias in Body axes [m/s^2]
            bias = obj.acc_Bias;
        end
        
        function setAccEstBias(obj, bias) 
            %% SETACCESTBIAS: sets accelerometer non-EKF estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % inputs: bias - 3x1 Real - estimated Bias in Body axes [m/s^2]
            obj.acc_Bias = bias;
        end
        
        function addAccEstBias(obj, dbias) 
            %% ADDACCESTBIAS: adds value to accelerometer non-EKF estimated bias
            %
            % the total estimated bias is composed by two components: one
            % is hard-compensated in all sensor measurements (non-EKF), while 
            % the other is present inside the EKF estimated state. Sensor
            % calibration is performed by transporting the EKF estimate to
            % the non-EKF estimate. Notice that this has no effect in the
            % total bias estimate, but it has an effect on the non-EKF
            % compensated estimates given by the CF.
            %
            % inputs: incremental dbias - 3x1 Real - incremental Bias [m/s^2]
            obj.acc_Bias = obj.acc_Bias + dbias;
        end
        
        function setLLA0(obj, lat, lon, h)
            %% SETLLA0: sets LLA0 WGS-84 reference position
            %
            % these LLA coordinates are supposed to be set to the initial
            % drone position. After initial configuration, these
            % coordinates are not supposed to be changed or updated in case
            % the drone moves. instead, the NED position vector is used to
            % define the relative position between LLA0 and the actual
            % drone position. here, obviously, we assume small planet
            % displacements.
            %
            % inputs: lat - Real scalar - latitude [rad]
            %         lon - Real scalar - longitude [rad]
            %         h - Real scalar - altitude [m]
            obj.LLA0 = [lat; lon; h];
        end
        
        function [lat, lon, h] = getLLA0(obj)
            %% GETLLA0: gets LLA0 WGS-84 reference position
            % 
            % these LLA coordinates are supposed to be set to the initial
            % drone position. After initial configuration, these
            % coordinates are not supposed to be changed or updated in case
            % the drone moves. instead, the NED position vector is used to
            % define the relative position between LLA0 and the actual
            % drone position. here, obviously, we assume small planet
            % displacements.
            %
            % outputs: lat - Real scalar - latitude [rad]
            %          lon - Real scalar - longitude [rad]
            %          h - Real scalar - altitude [m]
            lat = obj.LLA0(1);
            lon = obj.LLA0(2); 
            h = obj.LLA0(3);
        end
        
        function gain = getMagGainCF(obj)
            %% GETMAGGAINCF: gets complementary filter magnetic gain
            %
            % outputs: gain - Real scalar - CF Magnetic gain
            gain = obj.CF_Km;
        end
        
        function setMagGainCF(obj, gain)
            %% SETMAGGAINCF: sets complementary filter magnetic gain
            %
            % inputs: gain - Real scalar - CF Magnetic gain
            obj.CF_Km = gain;
        end
        
        function gain = getGravGainCF(obj)
            %% GETGRAVGAINCF: gets complementary filter gravitational gain
            %
            % outputs: gain - Real scalar - CF Gravitational gain
            gain = obj.CF_Kg;
        end
        
        function setGravGainCF(obj, gain)
            %% SETGRAVGAINCF: sets complementary filter gravitational gain
            %
            % inputs: gain - Real scalar - CF Gravitational gain
            obj.CF_Kg = gain;
        end
        
        function obj=lcNavi()
            %% lcNavi: class constructor
            
            % Empty constructor
        end
        
        function setTempo(obj, t)
            %% SETTEMPO: sets tempo in secs
            %
            % every update navigation routine in this library that asks for
            % the time uses the given data and automatically updates the
            % current time. therefore, getTempo() and setTempo() should not
            % be used by end-user. this function is called on
            % update/predict functions only. user should not care about
            % this.
            %
            %  inputs: t - Scalar Real - current time [sec]
            
            obj.tempo = t;
        end
        
    end
    
    %% Navigator static methods
    methods (Access=private)
        
        function Hk = getGNSSlooselyPosModel(obj)
            %% GETGNSSLOOSELYPOSMODEL: constructs Hk matrix for loosely coupled
            %
            % outputs: Hk - 3x18 Real - GNSS loosely coupled sensor model
            Hk = zeros(3,18);
            Hk(1:3,1:3) = eye(3);
        end
        
        function Hk = getGNSSlooselyVelModel(obj)
            %% GETGNSSLOOSELYVELMODEL: constructs Hk matrix for loosely coupled
            %
            % outputs: Hk - 3x18 Real - GNSS loosely coupled sensor model
            Hk = zeros(3,18);
            Hk(1:3,4:6) = eye(3);
        end
        
        function [Fk, Gk] = getINSErrorModel(obj, gyr, acc, mag, dt)
            %% GETINSERRORMODEL: constructs Fk matrix from INS error model
            %
            % inputs: gyr - 3x1 Real - uncompensated gyro measurement [rad/s]
            %         acc - 3x1 Real - uncompensated accelerometer measurement [m/s^2]
            %         mag - 3x1 Real - uncompensated magnetometer measurement [T]
            %         dt - Scalar Real - time difference since last update [sec]
            % outputs: Fk - 18x18 Real - discrete process matrix
            %          Gk - 18x12 Real - discrete process matrix
            
            %% matrices initialization
            At = zeros(18,18);
            Gk = zeros(18,12);
            %% position error model
            At(1:3,4:6) = eye(3);
            %% velocity error model
            D = obj.getAttitude(obj.MODE_DCM);
            accf = acc - obj.getAccEstBias();
            At(4:6,7:9) = D'*obj.skew(accf);
            At(4:6,10:12) = D';
            %% "psi-angle" error model
            wf = obj.angVelCF(gyr, acc, mag);
            km = obj.getMagGainCF();
            ka = obj.getGravGainCF();
            [inc, dec, amp] = obj.getLocalMag();
            Bi = obj.magVecNED(inc, dec, amp);
            Bc = D*Bi;
            gc = D*[0;0;obj.G]; 
            At(7:9,7:9) = -obj.skew(wf) + km*(obj.skew(Bc))^2 + ka*(obj.skew(gc))^2;
            At(7:9,10:12) = -ka*obj.skew(gc);
            At(7:9,13:15) = -eye(3);
            At(7:9,16:18) = km*obj.skew(Bc);
            %% velocity noise model
            Gk(4:6,1:3) = D';
            %% "psi-angle" noise model
            Gk(7:9,1:3) = -ka*obj.skew(gc);
            Gk(7:9,4:6) = -eye(3);
            Gk(7:9,7:9) = km*obj.skew(Bc);
            Gk(7:9,10:12) = -ka*obj.skew(gc);
            %% Finally, we discretize continuous error model
            Ts = dt;
            Fk = eye(18) + At*Ts;
        end
        
        function updateCF(obj, gyr, acc, mag, dt)
            %% UPDATECF: updates the estimate of the complementary filter (CF)
            %
            % inputs: gyr - 3x1 Real - uncompensated gyro measurement [rad/s]
            %         acc - 3x1 Real - uncompensated accelerometer measurement [m/s^2]
            %         mag - 3x1 Real - uncompensated magnetometer measurement [T]
            %         dt - Real scalar - time difference since last update [sec]
        
            % computes CF angular velocity
            wf = obj.angVelCF(gyr, acc, mag);
            % computes quaternion time derivative
            Qa = [0 -wf'; wf -obj.skew(wf)];
            q = obj.getAttitude(obj.MODE_QUATERNION);
            dq_dt = 1/2*Qa*q;
            % computes navigation derivatives
            D = obj.getAttitude(obj.MODE_DCM);
            accf = acc - obj.getAccEstBias();
            dv_dt = D'*accf+[0;0;obj.G];
            dp_dt = obj.getVelocity();
            % update navigation variables
            Ts = dt; % sampling time querying
            obj.setAttitude(obj.MODE_QUATERNION,quatnormalize((q+dq_dt*Ts)')');
            obj.setVelocity(dp_dt + dv_dt*Ts);
            pos = obj.getPosition();
            obj.setPosition(pos + dp_dt*Ts);
        end
        
        function wf = angVelCF(obj, gyr, acc, mag)
            %% ANGVELCF: gets complementary filter (CF) angular velocity
            %
            % inputs: gyr - 3x1 Real - uncompensated gyro measurement [rad/s]
            %         acc - 3x1 Real - uncompensated accelerometer measurement [m/s^2]
            %         mag - 3x1 Real - uncompensated magnetometer measurement [T]
            % outputs: wf - 3x1 Real - CF angular velocity [rad/s]
        
            %% compensates for bias in magnetometer
            magf = mag - obj.getMagEstBias();
            %% we are only interested in magnetic direction
            if (norm(magf,2)>0)
                B_body = magf;
            else
                error('Something is stinky at magnetometer measurements: Zero field magnitude');
            end
            %% get DCM attitude from last step
            D = obj.getAttitude(obj.MODE_DCM);
            %% compare gravity to accelerometer measurement
            % compensate for accelerometer bias
            accf = acc - obj.getAccEstBias();
            % I know, I know, here we assume equilibrium flight;
            % The beauty of this algorithm is that the Kalman filter will
            % also know, so it is not a big deal and sensors are low-cost
            % anyway!
            nG3 = cross(-accf,D*[0;0;obj.G]);
            %% compare magnetic field to magnetometer measurement
            [inc, dec, amp] = obj.getLocalMag();
            B_dir = obj.magVecNED(inc, dec, amp);
            nM3 = cross(B_body, D*B_dir);
            %% compute filtered angular velocity of body wrt inertial
            % compensate for gyro bias
            gyrf = gyr - obj.getGyroEstBias();
            % compute complementary filter
            wf = gyrf + obj.getGravGainCF()*nG3 + obj.getMagGainCF()*nM3;
        end
        
        function M = skew(obj,w)
            %% SKEW: returns matrix operator for the vector product with w
            %
            % input: w - 3x1 Real
            % output: M - 3x3 Real
            M = [ 0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0 ];
        end
        
        function B = magVecNED(obj,inclination, declination, magnit)
            %% MAGVECNED: returns magnetic vector in NED coordinates
            %
            % input: inclination - Real scalar [rad]
            %        declination - Real scalar [rad]
            %        magnit - Real scalar [Tesla]
            % output: B - 3x1 Real [Tesla]
            %
            % ATENTION: mag field pointing downwards has positive inclination!
            % while the sign of the declination is the same as defined in
            % [www.magnetic-declination.com].
            B = magnit*[cos(inclination)*cos(declination); ...
                cos(inclination)*sin(declination); ...
                sin(inclination) ];
        end
        
        function setOpMode(obj, mode)
            %% SETOPMODE: sets INS operation mode
            %
            % inputs: mode - Real scalar - constant mode
            % NOTICE: the constants can be OP_IDLE, OP_STATIC_CALIBRATION,
            % OP_DYNAMIC_CALIBRATION, OP_CRUISE
            if ( mode == obj.OP_IDLE || mode == obj.OP_STATIC_CALIBRATION || mode == obj.OP_DYNAMIC_CALIBRATION || mode == obj.OP_CRUISE || mode == obj.OP_ALIGNMENT )
                obj.OpMode = mode;
                return
            end
            error('Invalid INS operation mode!');
        end
        
        function mode = getOpMode(obj)
            %% GETOPMODE: gets INS operation mode
            %
            % outputs: mode - Real scalar - constant mode
            % NOTICE: the constants can be OP_IDLE, OP_MAG_CALIBRATION,
            % OP_GYR_CALIBRATION, OP_CRUISE
            mode = obj.OpMode;
        end
        
        function pos_ned = lla2posNED(obj, lla)
            %% GETLLA: gets NED coordinates wrt LLA0 of given arbitrary LLA
            %
            % this function performs the computation to compute reference
            % differential NED position from arbitrary given LLA and drone LLAO reference
            % position. this function is necessary because we keep position
            % information as NED position vector and eventually LLA
            % coordinates are required instead. NOTICE: the LLAO component
            % in this class is a reference point for the zero-position of
            % the NED position cordinates, and NOT the LLA coordinate of
            % the drone.
            %
            % outputs : pos_ned - 3x1 Real - arbitrary NED position coordinates [m,m,m]
            
            %% first conversion of WGS-84 to local NED
            lat = lla(1);
            lon = lla(2);
            h = lla(3);
            [lat0, lon0, h0] = obj.getLLA0();
            %% get differentials in LLA
            dlat = lat - lat0;
            dlon = lon - lon0;
            dh   = h - h0; 
            %% transform the LLA differentials in NED differentials
            a = 6378137; % semi-major axis of Earth (in meters) 
            ec2 = 6.69437999014e-3; % first eccentricity squared
            ec = sqrt(ec2); % first eccentricity
            RN = a*(1-ec^2)/sqrt((1-(ec*sin(lat0))^2)^3); 
            RE = a/sqrt(1-(ec*sin(lat0))^2);
            pn_gnss = dlat*(RN + h0);
            pe_gnss = dlon*(RE + h0)*cos(lat0);
            pd_gnss = -dh;
            pos_ned = [pn_gnss; pe_gnss; pd_gnss];
            
        end
        
        function t = getTempo(obj)
            %% GETTEMPO: gets tempo in secs
            %
            % every update navigation routine in this library that asks for
            % the time uses the given data and automatically updates the
            % current time. therefore, getTempo() and setTempo() should not
            % be used by end-user. this function is called on
            % update/predict functions only. user should not care about
            % this.
            %
            %  outputs: t - Scalar Real - current time [sec]
            t = obj.tempo;
        end
        
    end
    
end



