﻿#ifndef CUSTOMMATH_H
#define CUSTOMMATH_H

#define DEG2RAD (0.01745329251994329576923690768489)
#define RAD2DEG 1/DEG2RAD
#define GRAVITY 9.80665
#define PI 3.1415926535897932384626433

#include "math.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <string.h>
using namespace Eigen;

#define CustomArraySize 50
namespace Eigen
{
		typedef Matrix<double, -1, -1, 0, CustomArraySize, CustomArraySize> MatrixCXd;
		typedef Matrix<double, -1, -1, 0, CustomArraySize, 1> VectorCXd;
        typedef Matrix<int, -1, -1, 0, CustomArraySize, CustomArraySize> MatrixCXi;
		typedef Matrix<int, -1, -1, 0, CustomArraySize, 1> VectorCXi;
        typedef Matrix<bool, -1, -1, 0, CustomArraySize, CustomArraySize> MatrixCXb;
		typedef Matrix<bool, -1, -1, 0, CustomArraySize, 1> VectorCXb;
}

namespace CustomMath
{

    // pseudo inverse
    static MatrixXd pseudoInverseSVD(const MatrixXd& a, double epsilon = std::numeric_limits<double>::epsilon())
    {
        JacobiSVD<MatrixXd> svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
        double tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);

        MatrixXd ainv(a.cols(), a.rows());
        ainv.noalias() = svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();

        return ainv;
    }

    static MatrixXd pseudoInverseQR(const MatrixXd& A)
    {
        CompleteOrthogonalDecomposition<MatrixXd> cod(A);
        return cod.pseudoInverse();
    }

    static MatrixXd OneSidedInverse(const MatrixXd& A) //no use
    {
        MatrixXd Ainv(A.cols(), A.rows());

        if (A.cols() > A.rows()) //left inverse
        {
            MatrixXd AAT(A.rows(), A.rows());
            AAT = A * A.transpose();
            MatrixXd AATinv(A.rows(), A.rows());
            AATinv = AAT.inverse();

            Ainv = A.transpose() * AATinv;
        }
        else //right inverse
        {
            MatrixXd ATA(A.cols(), A.cols());
            ATA = A.transpose() * A;
            MatrixXd ATAinv(A.cols(), A.cols());
            ATAinv = ATA.inverse();
            Ainv = ATAinv * A.transpose();
        }

        return Ainv;
    }

    static MatrixXd WeightedPseudoInverse(const MatrixXd& Q, const MatrixXd& W, const bool Is_W_fullrank)
    {
        MatrixXd invQ(Q.cols(), Q.rows());
        MatrixXd transposeQ(Q.cols(), Q.rows());
        MatrixXd Winv(W.cols(), W.rows());
        if (Is_W_fullrank == true)
        {
            Winv = W.inverse();
        }
        else
        {
            Winv = pseudoInverseQR(W);
            //Winv = pseudoInverseSVD(W, 0.0001);
        }
        transposeQ = Q.transpose();
        MatrixXd tmp1(W.cols(), Q.rows());
        tmp1.noalias() = Winv * transposeQ;
        MatrixXd tmp2(Q.rows(), W.rows());
        tmp2.noalias() = Q * Winv;
        MatrixXd tmp3(Q.rows(), Q.rows());
        tmp3.noalias() = tmp2 * transposeQ;
        MatrixXd tmp3inv(Q.rows(), Q.rows());
        tmp3inv = tmp3.inverse();
        invQ.noalias() = tmp1 * tmp3inv;
        //invQ = Winv*transposeQ*pseudoInverseQR(Q*Winv*transposeQ);
        //invQ = Winv*transposeQ*pseudoInverseSVD(Q*Winv*transposeQ);

        return invQ;
    }

    static MatrixXd DampedWeightedPseudoInverse(const MatrixXd& Q, const MatrixXd& W, const bool Is_W_fullrank)
    {
        MatrixXd invQ(Q.cols(), Q.rows());
        MatrixXd transposeQ(Q.cols(), Q.rows());
        MatrixXd Winv(W.cols(), W.rows());
        if (Is_W_fullrank == true)
        {
            Winv = W.inverse();
        }
        else
        {
            Winv = pseudoInverseQR(W);
            //Winv = pseudoInverseSVD(W, 0.0001);
        }
        transposeQ = Q.transpose();
        double sigma = 0.01;
        MatrixXd Id(Q.rows(), Q.rows());
        Id.setIdentity();
        //invQ = Winv*transposeQ*pseudoInverseQR(Q*Winv*transposeQ + sigma*Id);
        MatrixXd tmp1(W.cols(), Q.rows());
        tmp1.noalias() = Winv * transposeQ;
        MatrixXd tmp2(Q.rows(), W.rows());
        tmp2.noalias() = Q * Winv;
        MatrixXd tmp3(Q.rows(), Q.rows());
        tmp3.noalias() = tmp2 * transposeQ;
        MatrixXd tmp4(Q.rows(), Q.rows());
        tmp4.noalias() = tmp3 + sigma * Id;
        MatrixXd tmp4inv(Q.rows(), Q.rows());
        //tmp4inv = pseudoInverseSVD(tmp4, 0.0001);
        tmp4inv = pseudoInverseQR(tmp4);

        invQ.noalias() = tmp1 * tmp4inv;

        return invQ;
    }

    // lowpass filter
    static double VelLowpassFilter(double rT, //delT
        double rWn, //cutoff freq (rad/sec)
        double rQ_pre, //previous position
        double rQ, //current position
        double rVel_pre //previous filtered velocity
    )
    {
        double rA;
        double rB;
        double rC;
        double rD;

        double rVel;

        rA = 2.0 * rWn;
        rB = (-2.0) * rWn;
        rC = 2.0 + rT * rWn;
        rD = rT * rWn - 2.0;

        rVel = ((-(rD)) * rVel_pre + (rA)*rQ + (rB)*rQ_pre) / (rC);

        return(rVel);
    }

    static double LowPassFilter(double dT, double Wc, double X, double preY) //sampling time, cutoff freq, input, previous output
    {
        double tau = 1.0 / Wc;
        double y = tau / (tau + dT) * preY + dT / (tau + dT) * X;
        return y;
    }



    static Eigen::Vector3d GetBodyRotationAngle(Eigen::Matrix3d RotMat)
    {
        double rollangle, pitchangle, yawangle;
        Eigen::Vector3d BodyAngle;
        BodyAngle.setZero();

        //  pitchangle = atan2(-RotMat(2,0),sqrt(RotMat(0,0)*RotMat(0,0)+RotMat(1,0)*RotMat(1,0)));
        //  yawangle = atan2(RotMat(1,0)/cos(pitchangle),RotMat(0,0)/cos(pitchangle));
        //  rollangle = atan2(RotMat(2,1)/cos(pitchangle), RotMat(2,2)/cos(pitchangle));

        double threshold = 0.001;
        pitchangle = -asin(RotMat(2, 0));
        if (RotMat(2, 0) > 1.0 - threshold && RotMat(2, 0) < 1.0 + threshold) //when RotMat(2,0) == 1 == sin(pitch)=1 : pitch angle = 90 deg
        {//Gimbal lock, pitch = -90deg == cos(pitch) = 0 : singularity is occured
            rollangle = atan2(-RotMat(0, 1), -RotMat(0, 2));
            yawangle = 0.0;
        }
        else if (RotMat(2, 0) < -1.0 + threshold && RotMat(2, 0) > -1.0 - threshold) //when RotMat(2,0) == -1 == sin(pitch)= -1 : pitch angle = -90 deg
        {//Gimbal lock, pitch = 90deg
            rollangle = atan2(RotMat(0, 1), RotMat(0, 2));
            yawangle = 0.0;
        }
        else //general solution
        {
            rollangle = atan2(RotMat(2, 1), RotMat(2, 2));
            yawangle = atan2(RotMat(1, 0), RotMat(0, 0));
        }

        BodyAngle(0) = rollangle;
        BodyAngle(1) = pitchangle;
        BodyAngle(2) = yawangle;

        return BodyAngle;
    }

    static double GetBodyPitchAngle(Eigen::Matrix3d RotMat)
    {
        double pitchangle;
        pitchangle = -asin(RotMat(2, 0));

        return pitchangle;
    }

    static double GetBodyRollAngle(Eigen::Matrix3d RotMat)
    {
        double rollangle, pitchangle;
        double threshold = 0.001;
        pitchangle = -asin(RotMat(2, 0));
        if (RotMat(2, 0) > 1.0 - threshold && RotMat(2, 0) < 1.0 + threshold) //when RotMat(2,0) == 1
        {//Gimbal lock, pitch = -90deg
            rollangle = atan2(-RotMat(0, 1), -RotMat(0, 2));
        }
        else if (RotMat(2, 0) < -1.0 + threshold && RotMat(2, 0) > -1.0 - threshold) //when RotMat(2,0) == -1
        {//Gimbal lock, pitch = 90deg
            rollangle = atan2(RotMat(0, 1), RotMat(0, 2));
        }
        else //general solution
        {
            rollangle = atan2(RotMat(2, 1), RotMat(2, 2));
        }

        return rollangle;
    }

    static double GetBodyYawAngle(Eigen::Matrix3d RotMat)
    {
        double pitchangle, yawangle;

        double threshold = 0.001;
        pitchangle = -asin(RotMat(2, 0));
        if (RotMat(2, 0) > 1.0 - threshold && RotMat(2, 0) < 1.0 + threshold) //when RotMat(2,0) == 1
        {//Gimbal lock, pitch = -90deg
            yawangle = 0.0;
        }
        else if (RotMat(2, 0) < -1.0 + threshold && RotMat(2, 0) > -1.0 - threshold) //when RotMat(2,0) == -1
        {//Gimbal lock, pitch = 90deg
            yawangle = 0.0;
        }
        else //general solution
        {
            yawangle = atan2(RotMat(1, 0), RotMat(0, 0));
        }

        return yawangle;
    }

    static Eigen::Matrix3d GetBodyRotationMatrix(double Roll, double Pitch, double Yaw)
    {
        Eigen::Matrix3d R_yaw;
        R_yaw.setZero();
        //yawŽÂ zÃà¿¡ ŽëÇÑ ÈžÀü
        R_yaw(2, 2) = 1.0;
        R_yaw(0, 0) = cos(Yaw);
        R_yaw(0, 1) = -sin(Yaw);
        R_yaw(1, 0) = sin(Yaw);
        R_yaw(1, 1) = cos(Yaw);

        Eigen::Matrix3d R_pitch;
        R_pitch.setZero();
        //pitchŽÂ yÃà¿¡ ŽëÇÑ ÈžÀü
        R_pitch(1, 1) = 1.0;
        R_pitch(0, 0) = cos(Pitch);
        R_pitch(2, 2) = cos(Pitch);
        R_pitch(0, 2) = sin(Pitch);
        R_pitch(2, 0) = -sin(Pitch);

        Eigen::Matrix3d R_roll;
        R_roll.setZero();
        //rollÀº xÃà¿¡ ŽëÇÑ ÈžÀü
        R_roll(0, 0) = 1.0;
        R_roll(1, 1) = cos(Roll);
        R_roll(2, 2) = cos(Roll);
        R_roll(1, 2) = -sin(Roll);
        R_roll(2, 1) = sin(Roll);

        Eigen::Matrix3d RGyro;
        RGyro.noalias() = R_yaw * R_pitch * R_roll;
        //dMatrix RGyro_inv(3,3);
        //RGyro.inv(RGyro_inv);

        //return RGyro_inv;
        return RGyro;
    }

    static Eigen::Matrix3d rotateWithZ(double yaw_angle)
    {
        Eigen::Matrix3d rotate_wth_z(3, 3);

        rotate_wth_z(0, 0) = cos(yaw_angle);
        rotate_wth_z(1, 0) = sin(yaw_angle);
        rotate_wth_z(2, 0) = 0.0;

        rotate_wth_z(0, 1) = -sin(yaw_angle);
        rotate_wth_z(1, 1) = cos(yaw_angle);
        rotate_wth_z(2, 1) = 0.0;

        rotate_wth_z(0, 2) = 0.0;
        rotate_wth_z(1, 2) = 0.0;
        rotate_wth_z(2, 2) = 1.0;

        return rotate_wth_z;
    }

    static Eigen::Matrix3d rotateWithY(double pitch_angle)
    {
        Eigen::Matrix3d rotate_wth_y(3, 3);

        rotate_wth_y(0, 0) = cos(pitch_angle);
        rotate_wth_y(1, 0) = 0.0;
        rotate_wth_y(2, 0) = -sin(pitch_angle);

        rotate_wth_y(0, 1) = 0.0;
        rotate_wth_y(1, 1) = 1.0;
        rotate_wth_y(2, 1) = 0.0;

        rotate_wth_y(0, 2) = sin(pitch_angle);
        rotate_wth_y(1, 2) = 0.0;
        rotate_wth_y(2, 2) = cos(pitch_angle);

        return rotate_wth_y;
    }

    static Eigen::Matrix3d rotateWithX(double roll_angle)
    {
        Eigen::Matrix3d rotate_wth_x(3, 3);

        rotate_wth_x(0, 0) = 1.0;
        rotate_wth_x(1, 0) = 0.0;
        rotate_wth_x(2, 0) = 0.0;

        rotate_wth_x(0, 1) = 0.0;
        rotate_wth_x(1, 1) = cos(roll_angle);
        rotate_wth_x(2, 1) = sin(roll_angle);

        rotate_wth_x(0, 2) = 0.0;
        rotate_wth_x(1, 2) = -sin(roll_angle);
        rotate_wth_x(2, 2) = cos(roll_angle);

        return rotate_wth_x;
    }

    static Eigen::Matrix3d skew(Eigen::Vector3d src)
    {
        Eigen::Matrix3d skew;
        skew.setZero();
        skew(0, 1) = -src[2];
        skew(0, 2) = src[1];
        skew(1, 0) = src[2];
        skew(1, 2) = -src[0];
        skew(2, 0) = -src[1];
        skew(2, 1) = src[0];

        return skew;
    }
    /*
    static Eigen::Vector3d getPhi(Eigen::Matrix3d current_rotation, Eigen::Matrix3d desired_rotation)
    {
        Eigen::Vector3d phi;
        Eigen::Vector3d s[3], v[3], w[3];

        for (int i = 0; i < 3; i++)
        {
            v[i] = current_rotation.block<3, 1>(0, i);
            w[i] = desired_rotation.block<3, 1>(0, i);
            s[i] = v[i].cross(w[i]);
        }
        phi = s[0] + s[1] + s[2];
        phi = -0.5 * phi;

        return phi;
    }
    */
    static Eigen::Vector3d getPhi(Eigen::Matrix3d& RotationMtx, Eigen::Matrix3d& DesiredRotationMtx) //Orientation 구성
    {
        //Get SkewSymmetric
        Eigen::Matrix3d s1_skew;
        Eigen::Matrix3d s2_skew;
        Eigen::Matrix3d s3_skew;

        Eigen::Vector3d RotMtxcol1;
        Eigen::Vector3d RotMtxcol2;
        Eigen::Vector3d RotMtxcol3;
        for (int i = 0; i < 3; i++)
        {
            RotMtxcol1(i) = RotationMtx(i, 0);
            RotMtxcol2(i) = RotationMtx(i, 1);
            RotMtxcol3(i) = RotationMtx(i, 2);
        }

        s1_skew = skew(RotMtxcol1);
        s2_skew = skew(RotMtxcol2);
        s3_skew = skew(RotMtxcol3);
        /////////////////////////////////////////////////////////////

        Eigen::Vector3d s1d;
        Eigen::Vector3d s2d;
        Eigen::Vector3d s3d;
        for (int i = 0; i < 3; i++)
        {
            s1d(i) = DesiredRotationMtx(i, 0);
            s2d(i) = DesiredRotationMtx(i, 1);
            s3d(i) = DesiredRotationMtx(i, 2);
        }

        Eigen::Vector3d s1f;
        Eigen::Vector3d s2f;
        Eigen::Vector3d s3f;

        s1f = s1_skew * s1d;
        s2f = s2_skew * s2d;
        s3f = s3_skew * s3d;
        /////////////////////////////////////////////////////////////
        //phi.resize(3);
        Eigen::Vector3d phi;
        phi = (s1f + s2f + s3f) * (-1.0 / 2.0);
        return phi;
    }

    static Eigen::Vector3d OrientationVelocity(Eigen::Matrix3d Rot, Eigen::Matrix3d Rotdot) //rotation matrix, derivative of rotation matrix
    {
        Eigen::Matrix3d RdotRT;
        RdotRT = Rotdot * Rot.transpose();
        Eigen::Vector3d OriVel;
        OriVel(0) = RdotRT(2, 1);
        OriVel(1) = RdotRT(0, 2);
        OriVel(2) = RdotRT(1, 0);

        return OriVel;
    }


    static double Cubic(double rT, double rT_0, double rT_f, double rx_0, double rx_dot_0, double rx_f, double rx_dot_f)
    {
        double rx_q;

        if (rT < rT_0)
        {
            rx_q = rx_0;
        }
        else if (rT > rT_f)
        {
            rx_q = rx_f;
        }
        else {

            rx_q = rx_0 + rx_dot_0 * (rT - rT_0)
                + (3.0 * (rx_f - rx_0) / ((rT_f - rT_0) * (rT_f - rT_0)) - 2.0 * rx_dot_0 / (rT_f - rT_0) - rx_dot_f / (rT_f - rT_0)) * (rT - rT_0) * (rT - rT_0)
                + (-2.0 * (rx_f - rx_0) / ((rT_f - rT_0) * (rT_f - rT_0) * (rT_f - rT_0)) + (rx_dot_0 + rx_dot_f) / ((rT_f - rT_0) * (rT_f - rT_0))) * (rT - rT_0) * (rT - rT_0) * (rT - rT_0);
        }
        return (rx_q);
    }

    static double CubicDot(double rT, double rT_0, double rT_f, double rx_0, double rx_dot_0, double rx_f, double rx_dot_f)
    {
        double rx_q_dot;

        if (rT < rT_0)
        {
            rx_q_dot = rx_dot_0;
        }
        else if (rT > rT_f)
        {
            rx_q_dot = rx_dot_f;
        }
        else {
            rx_q_dot = rx_dot_0 + 2.0 * (3.0 * (rx_f - rx_0) / ((rT_f - rT_0) * (rT_f - rT_0)) - 2.0 * rx_dot_0 / (rT_f - rT_0) - rx_dot_f / (rT_f - rT_0)) * (rT - rT_0)
                + 3.0 * (-2.0 * (rx_f - rx_0) / ((rT_f - rT_0) * (rT_f - rT_0) * (rT_f - rT_0)) + (rx_dot_0 + rx_dot_f) / ((rT_f - rT_0) * (rT_f - rT_0))) * (rT - rT_0) * (rT - rT_0);
        }
        return (rx_q_dot);
    }
    static VectorXd LinearInterpolation(VectorXd init_pos, VectorXd goal_pos, double cur_time, double init_time, double goal_time)
    {
        VectorXd _slope;
        VectorXd _desired;
        _slope.setZero(init_pos.rows());
        _desired.setZero(init_pos.rows());
        _slope = ((goal_pos - init_pos) / (goal_time - init_time));
        _desired = _slope*(cur_time - init_time) + init_pos;

        return _desired;
    }
    static VectorXd CubicSplineInterpolation(VectorXd Pos_0, VectorXd Vel_0, VectorXd Pos_f, VectorXd Vel_f, double T, double T_0, double T_f, double type)
    {
        VectorXd xd(Pos_0.rows());
        VectorXd _b0, _b1, _b2, _b3;
        _b0.setZero(Pos_0.rows());
        _b1.setZero(Pos_0.rows());
        _b2.setZero(Pos_0.rows());
        _b3.setZero(Pos_0.rows());
        double _dt, _duration_time;
        double _motion_threshold = 0.001;

        _dt = T - T_0;
        _duration_time = T_f - T_0;
        _b0 = Pos_0;
        _b1 = Vel_0;
        _b2 = (3*(Pos_f - Pos_0) - _duration_time*(2*Vel_0 + Vel_f))/pow(_duration_time,2);
        _b3 = (-2*(Pos_f - Pos_0) + _duration_time*(Vel_0 + Vel_f))/pow(_duration_time,3);
        if(type == 0){ // pos
            if (T <= T_0)
            {
                xd = Pos_0;
            }
            else if (T >= T_f)
            {
                xd = Pos_f;
            }
            else{
                xd = _b0 + _b1*_dt + _b2*pow(_dt,2) + _b3*pow(_dt,3);
            }
            for (int i = 0; i < Pos_0.rows(); i++) //do not use cubic spline when desired motion is small
            {
                if (abs(Pos_f(i) - Pos_0(i)) <= _motion_threshold)
                {
                    xd(i) = Pos_f(i);
                }
            }
        }
        else{//vel
            if (T <= T_0)
            {
                xd = Vel_0;
            }
            else if (T >= T_f)
            {
                xd = Vel_f;
            }
            else{
                xd = _b1 + 2*_b2*_dt + 3*_b3*pow(_dt,2);
            }
            for (int i = 0; i < Vel_0.rows(); i++) //do not use cubic spline when desired motion is small
            {
                if (abs(Pos_f(i) - Pos_0(i)) <= _motion_threshold)
                {
                    xd(i) = 0.0;
                }
            }
        }
        return xd;
    }

    static double Min(double val1, double val2)
    {
        double min_val = 0.0;

        if (val1 >= val2)
        {
            min_val = val2;
        }
        else
        {
            min_val = val1;
        }

        return min_val;
    }

    static double Max(double val1, double val2)
    {
        double min_val = 0.0;

        if (val1 >= val2)
        {
            min_val = val1;
        }
        else
        {
            min_val = val2;
        }

        return min_val;
    }

    static double SwitchFunction(double start, double end, double val)
    {
        double alpha = 0.0;
        alpha = Cubic(val, start, end, 0.0, 0.0, 1.0, 0.0);
        return alpha;
    }

    static double linear(double x_0, double x_f, double xdot_0, double xdot_f, double t_0, double t_f, double t)
    {
        double slope = (x_f - x_0)/(t_f - t_0);
        double desired = slope*(t - t_0) + x_0;
        return desired;
    }
    static Matrix3d rotation(Vector3d x_0, Vector3d x_f, double t, double t_0, double t_f)
    {
        Matrix3d rotation_0, rotation_f;
        rotation_0 = GetBodyRotationMatrix(x_0(0), x_0(1), x_0(2));
        rotation_f = GetBodyRotationMatrix(x_f(0), x_f(1), x_f(2));
        double tau = linear(0,1,0,0, t_0, t_f, t);

        Matrix3d rot_scaler_skew;
        rot_scaler_skew = (rotation_0.transpose() * rotation_f).log();
        Matrix3d result = rotation_0 * (rot_scaler_skew * tau).exp();
        return result;
        // return rot_scaler_skew;
    }
    static Vector3d rotationDot(Vector3d x_0, Vector3d x_f, double t, double t_0, double t_f)
    {
        Matrix3d rotation_0, rotation_f;
        rotation_0 = GetBodyRotationMatrix(x_0(0),x_0(1),x_0(2));
        rotation_f = GetBodyRotationMatrix(x_f(0),x_f(1),x_f(2));

        Matrix3d rot_scaler_skew;
        rot_scaler_skew = (rotation_0.transpose() * rotation_f).log();
        Vector3d a, b, c, r;
	    double tau = (t - t_0) / (t_f - t_0);
        r(0) = rot_scaler_skew(2, 1);
        r(1) = rot_scaler_skew(0, 2);
        r(2) = rot_scaler_skew(1, 0);
        Vector3d rd;
        for (int i = 0; i < 3; i++)
        {
            rd(i) = linear(0, r(i), 0, 0, t_0, t_f, t);
        }
        rd = rotation_0 * rd;

        return rd;
    }



    static void ContactWrenchConeConstraintGenerate(double friction_coeff, double cop_x_b, double cop_y_b, double pushing_force_threshold, Eigen::MatrixXd& A, Eigen::VectorXd& lb, Eigen::VectorXd& ub) //A: constraint matrix
    {
        //This function is for contact wrench cone constraint
        //based on Caron, Stéphane, Quang-Cuong Pham, and Yoshihiko Nakamura. "Stability of surface contacts for humanoid robots: Closed-form formulae of the contact wrench cone for rectangular support areas." 2015 IEEE International Conference on Robotics and Automation (ICRA). IEEE, 2015.
        // * z direction denotes normal direction of contact plane (negative direction = toward object or ground)
        // * pushing_force_thrreshold should be negative number 
        // * coordinate has to locate center of the contact plane (cop_x_b = l_x/2, cop_y_b = l_y/2)
        // * friction_coeff, cop_x_b and cop_y_b should be positive number
        
        A.setZero(17, 6);
        lb.setZero(17);
        ub.setZero(17);
        //  A.setZero(16, 6);
        // lb.setZero(16);
        // ub.setZero(16);

        double inf_val = 100000000.0; //large number to describe infinite

        // A matrix
        // //CoP x
        // A(0, 2) = -cop_x_b;
        // A(0, 4) = -1.0;
        // A(1, 2) = cop_x_b;
        // A(1, 4) = -1.0;
        // //CoP y
        // A(2, 2) = -cop_y_b;
        // A(2, 3) = 1.0;
        // A(3, 2) = cop_y_b;
        // A(3, 3) = 1.0;
        // //Fx friction
        // A(4, 0) = 1.0;
        // A(4, 2) = -friction_coeff;
        // A(5, 0) = 1.0;
        // A(5, 2) = friction_coeff;
        // //Fy friction
        // A(6, 1) = 1.0;
        // A(6, 2) = -friction_coeff;
        // A(7, 1) = 1.0;
        // A(7, 2) = friction_coeff;
        // //Mz friction
        // A(8, 0) = -cop_y_b;
        // A(8, 1) = -cop_x_b;
        // A(8, 2) = -friction_coeff*(cop_x_b+cop_y_b)/2.0;
        // A(8, 3) = friction_coeff;
        // A(8, 4) = friction_coeff;
        // A(8, 5) = 1.0;
        // A(9, 0) = -cop_y_b;
        // A(9, 1) = cop_x_b;
        // A(9, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        // A(9, 3) = friction_coeff;
        // A(9, 4) = -friction_coeff;
        // A(9, 5) = 1.0;
        // A(10, 0) = cop_y_b;
        // A(10, 1) = -cop_x_b;
        // A(10, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        // A(10, 3) = -friction_coeff;
        // A(10, 4) = friction_coeff;
        // A(10, 5) = 1.0;
        // A(11, 0) = cop_y_b;
        // A(11, 1) = cop_x_b;
        // A(11, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        // A(11, 3) = -friction_coeff;
        // A(11, 4) = -friction_coeff;
        // A(11, 5) = 1.0;
        // A(12, 0) = -cop_y_b;
        // A(12, 1) = -cop_x_b;
        // A(12, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        // A(12, 3) = -friction_coeff;
        // A(12, 4) = -friction_coeff;
        // A(12, 5) = -1.0;
        // A(13, 0) = -cop_y_b;
        // A(13, 1) = cop_x_b;
        // A(13, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        // A(13, 3) = -friction_coeff;
        // A(13, 4) = friction_coeff;
        // A(13, 5) = -1.0;
        // A(14, 0) = cop_y_b;
        // A(14, 1) = -cop_x_b;
        // A(14, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        // A(14, 3) = friction_coeff;
        // A(14, 4) = -friction_coeff;
        // A(14, 5) = -1.0;
        // A(15, 0) = cop_y_b;
        // A(15, 1) = cop_x_b;
        // A(15, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        // A(15, 3) = friction_coeff;
        // A(15, 4) = friction_coeff;
        // A(15, 5) = -1.0;
        // //Fz (pushing force)
        // A(16, 2) = 1.0;
        
        // // lb vector        
        // lb(1) = -inf_val;
        // lb(3) = -inf_val;
        // lb(5) = -inf_val;
        // lb(7) = -inf_val;
        // lb(16) = -100; 
        
        // //ub vector
        // ub(0) = inf_val;
        // ub(2) = inf_val;
        // ub(4) = inf_val;
        // ub(6) = inf_val;
        // ub(8) = inf_val;
        // ub(9) = inf_val;
        // ub(10) = inf_val;
        // ub(11) = inf_val;
        // ub(12) = inf_val;
        // ub(13) = inf_val;
        // ub(14) = inf_val;
        // ub(15) = inf_val;
        // ub(16) = pushing_force_threshold; 

        
        // Fx friction
        A(0, 0) = -1.0;
        A(0, 2) = -friction_coeff;
        A(1, 0) = 1.0;
        A(1, 2) = -friction_coeff;
        //Fy friction
        A(2, 1) = -1.0;
        A(2, 2) = -friction_coeff;
        A(3, 1) = 1.0;
        A(3, 2) = -friction_coeff;
        //CoP x 
        A(4, 2) = -cop_y_b;
        A(4, 3) = -1.0;
        A(5, 2) = -cop_y_b;
        A(5, 3) = 1.0;
        //CoP y 
        A(6, 2) = -cop_x_b;
        A(6, 4) = -1.0;
        A(7, 2) = -cop_x_b;
        A(7, 4) = 1.0;
        //Mz friction
        A(8, 0) = -cop_y_b;
        A(8, 1) = -cop_x_b;
        A(8, 2) = -friction_coeff*(cop_x_b+cop_y_b);
        A(8, 3) = friction_coeff;
        A(8, 4) = friction_coeff;
        A(8, 5) = -1.0;
        A(9, 0) = -cop_y_b;
        A(9, 1) = cop_x_b;
        A(9, 2) = -friction_coeff * (cop_x_b + cop_y_b);
        A(9, 3) = friction_coeff;
        A(9, 4) = -friction_coeff;
        A(9, 5) = -1.0;
        A(10, 0) = cop_y_b;
        A(10, 1) = -cop_x_b;
        A(10, 2) = -friction_coeff * (cop_x_b + cop_y_b);
        A(10, 3) = -friction_coeff;
        A(10, 4) = friction_coeff;
        A(10, 5) = -1.0;
        A(11, 0) = cop_y_b;
        A(11, 1) = cop_x_b;
        A(11, 2) = -friction_coeff * (cop_x_b + cop_y_b);
        A(11, 3) = -friction_coeff;
        A(11, 4) = -friction_coeff;
        A(11, 5) = -1.0;
        A(12, 0) = cop_y_b;
        A(12, 1) = cop_x_b;
        A(12, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        A(12, 3) = friction_coeff;
        A(12, 4) = friction_coeff;
        A(12, 5) = 1.0;
        A(13, 0) = cop_y_b;
        A(13, 1) = -cop_x_b;
        A(13, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        A(13, 3) = friction_coeff;
        A(13, 4) = -friction_coeff;
        A(13, 5) = 1.0;
        A(14, 0) = -cop_y_b;
        A(14, 1) = cop_x_b;
        A(14, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        A(14, 3) = -friction_coeff;
        A(14, 4) = friction_coeff;
        A(14, 5) = 1.0;
        A(15, 0) = -cop_y_b;
        A(15, 1) = -cop_x_b;
        A(15, 2) = -friction_coeff * (cop_x_b + cop_y_b) / 2.0;
        A(15, 3) = -friction_coeff;
        A(15, 4) = -friction_coeff;
        A(15, 5) = 1.0;
        //Fz (pushing force)
        A(16, 2) = 1.0;
        if(pushing_force_threshold>0)
        {
        // // lb vector
        // lb(4) = -inf_val;  
        // lb(5) = -inf_val;  
        // lb(6) = -inf_val;  
        // lb(7) = -inf_val;        
        // lb(16) = -inf_val;
        
        // //ub vector
        // ub(0) = inf_val;
        // ub(1) = inf_val;
        // ub(2) = inf_val;
        // ub(3) = inf_val;
        // ub(8) = inf_val;
        // ub(9) = inf_val;
        // ub(10) = inf_val;
        // ub(11) = inf_val;
        // ub(12) = inf_val;
        // ub(13) = inf_val;
        // ub(14) = inf_val;
        // ub(15) = inf_val;
        // ub(16) = pushing_force_threshold;
            //lb vector
            lb(0) = -inf_val;
            lb(1) = -inf_val;
            lb(2) = -inf_val;
            lb(3) = -inf_val;
            lb(4) = -inf_val;  
            lb(5) = -inf_val;  
            lb(6) = -inf_val;  
            lb(7) = -inf_val; 
            lb(8) = -inf_val;
            lb(9) = -inf_val;
            lb(10) = -inf_val;
            lb(11) = -inf_val;
            lb(12) = -inf_val;
            lb(13) = -inf_val;
            lb(14) = -inf_val;
            lb(15) = -inf_val;       
            lb(16) = pushing_force_threshold;
            
            //ub vector
            ub(16) = inf_val;//pushing_force_threshold
        }
        else // right hand
        {
            //  lb vector   
            lb(8) = -inf_val;
            lb(9) = -inf_val;
            lb(10) = -inf_val;
            lb(11) = -inf_val;
            lb(12) = -inf_val;
            lb(13) = -inf_val;
            lb(14) = -inf_val;
            lb(15) = -inf_val;       
            lb(16) = -inf_val;
            
           // ub vector
            ub(0) = inf_val;
            ub(1) = inf_val;
            ub(2) = inf_val;
            ub(3) = inf_val; 
            ub(4) = inf_val;
            ub(5) = inf_val;  
            ub(6) = inf_val; 
            ub(7) = inf_val; 
            ub(16) = pushing_force_threshold;
        }
        // // lb vector
        // lb(4) = -inf_val;  
        // lb(5) = -inf_val;  
        // lb(6) = -inf_val;  
        // lb(7) = -inf_val;        
        // lb(16) = -inf_val;
        
        // //ub vector
        // ub(0) = inf_val;
        // ub(1) = inf_val;
        // ub(2) = inf_val;
        // ub(3) = inf_val;
        // ub(8) = inf_val;
        // ub(9) = inf_val;
        // ub(10) = inf_val;
        // ub(11) = inf_val;
        // ub(12) = inf_val;
        // ub(13) = inf_val;
        // ub(14) = inf_val;
        // ub(15) = inf_val;
        // ub(16) = pushing_force_threshold;
    }
    
    static void SwapFRandML(Eigen::VectorXd Force, VectorXd *SwapFRML)
    {   
        VectorXd _temp;
        _temp = Force.segment(3,3);
        Force.segment(3,3) = Force.segment(6,3);
        Force.segment(6,3) = _temp;

        *SwapFRML = Force;
        return;
    }
    

   
}





#endif

