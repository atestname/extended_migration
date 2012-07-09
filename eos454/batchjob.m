j = batch('cluster_part','configuration','torque2x4for24','matlabpool',7,'CaptureDiary',true)

% % I use 3 workers, but use your own
% % replace 'torque4x1' with the configuration you want
% % also in 'matlabpool' specify the max number of workers minus 1
% % (1 worker will be used as head node)
% 
% % IMPORTANT: remember your job id!
% j.state                         % shows what state the job is in
% j.display                       % shows job details
% j.id                            % this give the Job ID. Rememeber the Job ID!
% 
% % Quit and wait, and when coming back
% j = getClusterJob('torque4x1', 5)  % in this example my job ID is 5 (getClusterJob is a script I wrote)
% j.state                         % check that it says 'finished'
% data = j.getAllOutputArguments  % retrieve all the remote job variables into 'data'
% j.diary                         % displays all the command line window output of the remote job
% 
% % finally after you are all done:
% j.destroy                       % cleans up the job
