
function a = returnAction(p1)
    a = randsample([1, 2], 1, true, [p1 1-p1]);
end