function tests=test_deepsea3d_amplitude_classification
tests=functiontests(localfunctions);
end
function testUnstableReferenceNeverPhysicsFail(t)
m=struct('reflected_complex',.2,'increment_complex',.3);
s=classify_deepsea3d_amplitude(false,0,0,m);
verifyEqual(t,s,"INCONCLUSIVE");
end
function testNumericalFailureNeverPhysicsFail(t)
m=struct('reflected_complex',.2,'increment_complex',.3);
s=classify_deepsea3d_amplitude(true,.002,0,m);
verifyEqual(t,s,"INCONCLUSIVE");
end
function testIncrementCannotHideBehindBackground(t)
m=struct('reflected_complex',.0001,'increment_complex',.01);
s=classify_deepsea3d_amplitude(true,0,0,m);
verifyEqual(t,s,"CONDITIONAL_PHYSICAL_FAIL");
end
function testConditionalPass(t)
m=struct('reflected_complex',.0001,'increment_complex',.0002);
s=classify_deepsea3d_amplitude(true,0,0,m);
verifyEqual(t,s,"CONDITIONAL_G1_PASS");
end
