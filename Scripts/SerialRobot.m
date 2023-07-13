classdef SerialRobot

    properties
        DOF
        DH
        q
        T_E
        J_E
        centerOfMassRelativeToDH
        massOfTheLinks
        momentOfInertiaOfLinksTensor
        massMatrix
        gravitationalAcceleration

    end

    methods(Static)

        function translation_matrix = translation(L,axis)
            % TRANSLATION - Computes a 4x4 homogeneous translation matrix along a specified axis.
            %
            % Usage:
            %   translation_matrix = translation(L, axis)
            %
            % Inputs:
            %   L: Length or displacement along the specified axis
            %   axis: Axis of translation (1, 2, or 3) corresponding to x, y, or z-axis, respectively
            %
            % Outputs:
            %   translation_matrix: 4x4 symbolic matrix representing the translation
            %
            % Notes:
            %   - The function returns a symbolic matrix by default. If you require a
            %     numeric output, convert the output matrix to double using the 'double'
            %     function.
            %   - The axis parameter must be an integer between 1 and 3, corresponding to
            %     the x, y, or z-axis, respectively.
            %   - The translation matrix is a 4x4 homogeneous matrix with a one in the
            %     diagonal elements and zero elsewhere, except for the fourth column of
            %     the specified axis, which contains the translation distance L.
            %
            % Example:
            %   L = 2; % translation distance
            %   axis = 1; % x-axis
            %   T = translation(L, axis); % compute the translation matrix


            translation_matrix = sym(eye(4,4));

            translation_matrix(axis,4) = L;

        end

        function rotation_matrix = rotation( theta, axis)
            % ROTATION - Computes a 4x4 homogeneous rotation matrix for rotating about a coordinate axis.
            %
            % Usage:
            %   rotation_matrix = rotation(theta, axis)
            %
            % Inputs:
            %   theta: Angle of rotation in radians
            %   axis: Axis of rotation (1, 2, or 3) corresponding to x, y, or z-axis, respectively
            %
            % Outputs:
            %   rotation_matrix: 4x4 symbolic matrix representing the rotation
            %
            % Notes:
            %   - The function returns a symbolic matrix by default. If you require a
            %     numeric output, convert the output matrix to double using the 'double'
            %     function.
            %   - The axis parameter must be an integer between 1 and 3, corresponding to
            %     the x, y, or z-axis, respectively.
            %   - The rotation matrix is a 4x4 homogeneous matrix with a 3x3 rotation matrix
            %     in the upper-left corner, and a one in the diagonal elements and zero
            %     elsewhere, except for the specified axis column and row, which contains
            %     the rotation matrix elements.
            %
            % Example:
            %   theta = pi/4; % rotation angle
            %   axis = 2; % y-axis
            %   T = rotation(theta, axis); % compute the rotation matrix
            %
            % Author: [Your Name]
            % Contac

            rotation_matrix = sym( eye(4,4));

            if axis == 1

                R =[1         0              0
                    0,    cos(theta),  -sin(theta)
                    0,    sin(theta),   cos(theta)];
            elseif axis == 2
                R =[cos(theta)         0,        sin(theta)
                    0,              1,             0
                    -sin(theta),        0,        cos(theta)];

            elseif axis == 3
                R =[cos(theta),     -sin(theta),   0
                    sin(theta),     cos(theta),    0
                    0,              0 ,         1   ];

            end

            rotation_matrix(1:3,1:3)=R;

        end


        function T = transformMatrix(DH_param)
            % TRANSFORMMATRIX - Calculates the homogeneous transformation matrix from the end effector to the
            % inertia frame using Denavit-Hartenberg (DH) parameters in Craig's Robotic Book's convention.
            %
            % Usage:
            %   T = transformMatrix(DH_param)
            %
            % Inputs:
            %   DH_param: DH parameters matrix, where each row represents a joint. Each row should contain
            %             [alpha, a, d, theta] values for the joint in that order.
            %
            % Outputs:
            %   T: 4x4 symbolic matrix representing the transformation matrix
            %
            % Notes:
            %   - The function returns a symbolic matrix by default. If you require a
            %     numeric output, convert the output matrix to double using the 'double'
            %     function.
            %   - The DH_param matrix is a n x 4 matrix, where n is the number of joints in the
            %     robot. The four columns represent the DH parameters in the following order:
            %     alpha, a, d, and theta.
            %   - The transformation matrix is computed using the DH convention described in
            %     Craig's Robotic Book, which assumes that the coordinate frames are attached to
            %     the joints and the end effector frame is obtained by a sequence of rotations and
            %     translations.
            %
            % Example:
            %   DH_param = [alpha1, a1, d1, theta1;
            %               alpha2, a2, d2, theta2;
            %               ...
            %               alphan, an, dn, thetan]; % DH parameter matrix
            %
            %   T = transformMatrix(DH_param); % compute the transformation matrix
            transformation_matrix= sym(eye(4,4));

            DH_param_size = size(DH_param);
            DH_param_nrow = DH_param_size(1);

            for i = 1:DH_param_nrow

                transformation_matrix = transformation_matrix*SerialRobot.translation(DH_param(i,2),1)*SerialRobot.rotation(DH_param(i,1),1)*SerialRobot.rotation(DH_param(i,4),3)*SerialRobot.translation(DH_param(i,3),3);

            end

            T = simplify((transformation_matrix));

        end

        function inv_transform = inverseTrasformMatrix(tranformation_matrix)
            % INVERSETRANSFORMMATRIX - Computes the inverse of a given homogeneous
            % transformation matrix.
            %
            % Usage:
            %   inv_transform = inverseTransformMatrix(transformation_matrix)
            %
            % Inputs:
            %   transformation_matrix: 4x4 homogeneous transformation matrix to invert
            %
            % Outputs:
            %   inv_transform: 4x4 homogeneous transformation matrix representing the
            %   inverse transformation of the input matrix
            %
            % Notes:
            %   - The input transformation matrix must be a valid homogeneous
            %     transformation matrix, i.e., a 4x4 matrix with the bottom row consisting
            %     of [0 0 0 1].
            %   - The function returns a symbolic matrix by default. If you require a
            %     numeric output, convert the output matrix to double using the 'double'
            %     function.
            %
            % Example:
            %   T = [cos(theta) -sin(theta) 0 x;
            %        sin(theta) cos(theta)  0 y;
            %        0          0           1 z;
            %        0          0           0 1]; % homogeneous transformation matrix
            %
            %   inv_T = inverseTransformMatrix(T); % compute the inverse transformation
            %

            inv_transform = sym(eye(4,4));

            R_inv =  transpose(tranformation_matrix(1:3,1:3));
            inv_transform(1:3,1:3) = R_inv;

            inv_transform(1:3,4) =  - simplify(expand(R_inv* tranformation_matrix(1:3,4)));

        end

        function skewv = symskew(v)
            skewv = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ];
        end
        function exp_omega = angVelocity2Exp(omega, q)
            w = SerialRobot.symskew(omega);
            exp_omega = eye(3,3) + w*sin(q) + w*w*(1-cos(q));
            
        end
        function exp_coordinate = twist2Exp(screw, q)
            v = screw(1:3);
            w = screw(4:6);
            e_wt = SerialRobot.angVelocity2Exp(w, q);
            exp_coordinate = sym(eye(4,4));
            exp_coordinate(1:3,1:3) = e_wt;
            exp_coordinate(1:3,end) = (eye(3,3) - e_wt)*cross(w,v) ...
                + w*(dot(w, v))*q;
            if w == 0
                exp_coordinate(1:3,end) =v*q;
            end
        end
    end

    methods

        function obj = SerialRobot(q,DH_param,centerOfMassRelativeToDH,massOfTheLinks,momentOfInertiaOfLinksTensor,gravitationalAcceleration)
            % SERIALROBOT - A constructor function for a robotic manipulator object.
            %
            % Usage:
            %   obj = SerialRobot(q,DH_param,centerOfMassRelativeToDH,massOfTheLinks,momentOfInertiaOfLinksTensor)
            %
            % Inputs:
            %   q: A DOFx1 vector of current joint configurations of the robot, where DOF is the
            %      number of degrees of freedom of the robot.
            %
            %   DH_param: DH parameters matrix, where each row represents a joint. Each row should contain
            %             [alpha, a, d, theta] values for the joint in that order.
            %
            %   centerOfMassRelativeToDH (optional): A DOFx3 matrix representing the center of mass
            %                                        of each link with respect to their respective DH frames.
            %
            %   massOfTheLinks (optional): A DOFx1 vector representing the mass of each link.
            %
            %   momentOfInertiaOfLinksTensor (optional): A cell array of DOF matrices, where each matrix
            %                                            represents the moment of inertia tensor of each link
            %                                            with respect to their center of mass.
            %
            % Outputs:
            %   obj: A robotic manipulator object, which contains the following fields:
            %        - q: The current joint configurations of the robot.
            %        - DOF: The number of degrees of freedom of the robot.
            %        - DH: The DH parameters matrix of the robot.
            %        - centerOfMassRelativeToDH: The center of mass of each link with respect to their
            %                                    respective DH frames.
            %        - massOfTheLinks: The mass of each link.
            %        - momentOfInertiaOfLinksTensor: The moment of inertia tensor of each link with respect to
            %                                        their center of mass.
            %        - T_E: The homogeneous transformation matrix from the end effector to the
            %               robot base frame, calculated using the DH parameters.
            %        - J_E: The Jacobian matrix of the robot end effector, calculated using the DH parameters.
            %
            % Notes:
            %   - The function assumes that the robot's end effector is attached to the last joint of
            %     the robot and that the base frame is fixed.
            %   - If any of the optional input arguments are not provided, the corresponding fields in
            %     the object will be set to empty matrices.
            %
            % Example:
            %   DH_param = [alpha1, a1, d1, theta1;
            %               alpha2, a2, d2, theta2;
            %               ...
            %               alphan, an, dn, thetan]; % DH parameter matrix
            %   robot = SerialRobot(q,DH_param); % create a robot manipulator object
            %

            obj.q = q;

            obj.DOF = size(q,1);

            obj.DH = DH_param;

            if exist('centerOfMassRelativeToDH','var')

                obj.centerOfMassRelativeToDH = centerOfMassRelativeToDH;

            end

            if exist('massOfTheLinks','var')

                obj.massOfTheLinks = massOfTheLinks;

            end

            if exist('momentOfInertiaOfLinksTensor','var')

                obj.momentOfInertiaOfLinksTensor = momentOfInertiaOfLinksTensor;

            end

            if exist('gravitationalAcceleration','var')

                obj.gravitationalAcceleration = gravitationalAcceleration;
            else
                obj.gravitationalAcceleration = [0 0 -9.81]'
            end

            obj.T_E = SerialRobot.transformMatrix(DH_param);

            obj.J_E = obj.jacobianMatrix(DH_param);

        end


        function J = jacobianMatrix(obj,DH_parameter)
            % JACOBIANMATRIX - Calculates the Jacobian matrix for a robotic manipulator given its
            % Denavit-Hartenberg (DH) parameters and joint configurations.
            %
            % Usage:
            %   J = jacobianMatrix(obj,DH_parameter)
            %
            % Inputs:
            %   obj: An object of the robotic manipulator class, which contains the number of
            %        degrees of freedom (DOF) and the current joint configurations of the robot.
            %
            %   DH_parameter: DH parameters matrix, where each row represents a joint. Each row should contain
            %                 [alpha, a, d, theta] values for the joint in that order.
            %
            % Outputs:
            %   J: 6xDOF symbolic matrix representing the Jacobian matrix, where DOF is the number of
            %      degrees of freedom of the robot.
            %
            % Notes:
            %   - The function assumes that the robot's end effector is attached to the last joint of
            %     the robot and that the base frame is fixed.
            %   - The DH_parameter matrix is a n x 4 matrix, where n is the number of joints in the
            %     robot. The four columns represent the DH parameters in the following order:
            %     alpha, a, d, and theta.
            %   - The Jacobian matrix is computed using the chain rule and the partial derivatives
            %     of the transformation matrix with respect to each joint angle.
            %
            % Example:
            %   robot = RobotManipulator(DOF); % create a robot manipulator object
            %   DH_param = [alpha1, a1, d1, theta1;
            %               alpha2, a2, d2, theta2;
            %               ...
            %               alphan, an, dn, thetan]; % DH parameter matrix
            %   robot.q = [q1;q2;...;qDOF]; % set the joint configurations
            %   J = jacobianMatrix(robot,DH_param); % compute the Jacobian matrix

            transformation = obj.transformMatrix(DH_parameter);
            J_L = jacobian(transformation(1:3,end),obj.q);

            J_A = sym(zeros(3,obj.DOF));

            a = simplify(expand(transformation(2,1:3)));
            b = simplify(expand(jacobian(transformation(3,1:3),obj.q)));
            J_A(1,:) = a*b;

            a = simplify(expand(transformation(3,1:3)));
            b = simplify(expand(jacobian(transformation(1,1:3),obj.q)));
            J_A(2,:) = a*b;

            a = simplify(expand(transformation(1,1:3)));
            b = simplify(expand(jacobian(transformation(2,1:3),obj.q)));
            J_A(3,:) = a*b;

            J = simplify(expand([J_L;J_A]));

        end


        function M = massMatrixCalc(obj)
            % MASSMATRIXCALC - Calculates the mass matrix of a robotic manipulator given its
            % physical properties of the links.
            %
            % Usage:
            %   M = massMatrixCalc(obj)
            %
            % Inputs:
            %   obj: An object of the robotic manipulator class, which contains the number of
            %        degrees of freedom (DOF), the DH parameters, and the physical properties of the links.
            %
            % Outputs:
            %   M: A DOFxDOF symbolic matrix representing the mass matrix of the robot, where DOF
            %      is the number of degrees of freedom of the robot.
            %
            % Notes:
            %   - The function assumes that the robot's end effector is attached to the last joint of
            %     the robot and that the base frame is fixed.
            %   - The mass matrix is computed using the physical properties of the links, which include
            %     the center of mass of each link with respect to their respective DH frames, the mass
            %     of each link, and the moment of inertia tensor of each link with respect to their center
            %     of mass.
            %
            % Example:
            %   robot = SerialRobot(); % create a robot manipulator object
            %   robot.centerOfMassRelativeToDH = centerOfMassRelativeToDH; % set the center of mass of each link
            %   robot.massOfTheLinks = massOfTheLinks; % set the mass of each link
            %   robot.momentOfInertiaOfLinksTensor = momentOfInertiaOfLinksTensor; % set the moment of inertia tensor of each link
            %   M = massMatrixCalc(robot); % compute the mass matrix
            %
            M = sym(zeros(obj.DOF,obj.DOF));
            for i = 1:obj.DOF
                DH_param = [obj.DH(1:i,:); obj.centerOfMassRelativeToDH(i,:)];

                J_i = obj.jacobianMatrix(DH_param);
                J_Li = J_i(1:3,:);
                J_Ai = J_i(4:6,:);

                T_i = obj.transformMatrix(DH_param);
                R_i = T_i(1:3,1:3);
                I_tensor_i = R_i*obj.momentOfInertiaOfLinksTensor(:,:,i)*(R_i');

                m_i = obj.massOfTheLinks(i);

                M_i = simplify(m_i*(J_Li')*J_Li + (J_Ai')*I_tensor_i*J_Ai);
                M = M + M_i;

            end
            M = simplify(M);
        end


        function C = christoffelMatrixCalc(obj)

            DOF = obj.DOF;
            q = obj.q;
            C = sym(zeros(DOF,DOF,DOF));
            M = obj.massMatrixCalc();
            for i = 1:DOF
                for j = 1:DOF
                    for k = 1:DOF
                        C(i,j,k) = 1/2*(diff(M(k,j),q(i)) + diff(M(k,i),q(j))...
                            -diff(M(i,j),q(k)));
                    end
                end
            end

        end

        function C = velocityTermsDynamic(obj, q_dot)
            C = sym(zeros(obj.DOF,1));
            Christofel = obj.christoffelMatrixCalc();
            for i = 1:obj.DOF
                for j=1:obj.DOF
                    for k = 1:obj.DOF
                        C(i) = C(i) + Christofel(i,j,k)*q_dot(j)*q_dot(k);
                    end
                end
            end

        end
        function G = gravityEffectCalc(obj)
            % GRAVITYEFFECT - Calculates the gravity effect on a robotic manipulator given its
            % physical properties of the links and the acceleration due to gravity.
            %
            % Usage:
            %   G = gravityEffect(obj)
            %
            % Inputs:
            %   obj: An object of the robotic manipulator class, which contains the number of
            %        degrees of freedom (DOF), the DH parameters, and the physical properties of the links.
            %
            % Outputs:
            %   G: A DOFx1 symbolic vector representing the gravity effect on the robot, where DOF
            %      is the number of degrees of freedom of the robot.
            %
            % Notes:
            %   - The function assumes that the robot's end effector is attached to the last joint of
            %     the robot and that the base frame is fixed.
            %   - The gravity effect is computed using the physical properties of the links, which include
            %     the center of mass of each link with respect to their respective DH frames and the mass
            %     of each link, as well as the acceleration due to gravity.
            %   - The gravity effect is calculated using the Jacobian method, which involves calculating
            %     the geometric Jacobian of each link.
            %
            % Example:
            %   robot = SerialRobot(); % create a robot manipulator object
            %   robot.centerOfMassRelativeToDH = centerOfMassRelativeToDH; % set the center of mass of each link
            %   robot.massOfTheLinks = massOfTheLinks; % set the mass of each link
            %   G = gravityEffect(robot); % compute the gravity effect
            G = sym(zeros(obj.DOF,1));
            g = obj.gravitationalAcceleration;
            mJ = sym(zeros(3,obj.DOF));
            for i = 1:obj.DOF
                DH_param = [obj.DH(1:i,:); obj.centerOfMassRelativeToDH(i,:)];

                J_i = obj.jacobianMatrix(DH_param);
                J_Li = J_i(1:3,:);
                mJ = mJ + J_Li*obj.massOfTheLinks(i);
            end
            G = simplify(g'*mJ);

        end

        function q_ddot = forwardDynamics(obj)

            q_dot = sym('q_dot',[obj.DOF,1],'real');
            tau = sym('tau',[obj.DOF,1],'real');

            C = velocityTermsDynamic(obj,q_dot);

            q_ddot = massMatrixCalc(obj)\(tau - gravityEffectCalc(obj)- C);

            matlabFunction(q_ddot, 'File', 'forwardDynamics');


        end

        function tau = inverseDynamics(obj)
            q_ddot = sym('q_ddot', [obj.DOF,1], 'real');
            q_dot = sym('q_dot',[obj.DOF,1],'real');
            C = velocityTermsDynamic(obj,q_dot);

            tau = massMatrixCalc(obj)*q_ddot + C + gravityEffect(obj);

            matlabFunction(tau,'File','inverseDynamic')


        end
    end

end










