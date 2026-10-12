function evidence = fk_single_receiver_expect_error(fn,expected_id)
%FK_SINGLE_RECEIVER_EXPECT_ERROR An unrelated exception must fail validation.
actual='';
try
    fn();
catch exception
    actual=exception.identifier;
end
assert(strcmp(actual,expected_id),'fk_validation:UnexpectedError', ...
    'Expected error %s; actual error %s.',expected_id,actual);
evidence=struct('expected_id',expected_id,'actual_id',actual,'pass',true);
end
