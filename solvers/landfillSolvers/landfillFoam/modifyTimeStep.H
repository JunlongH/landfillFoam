// iter > 7 timestep = timestep * 0.7
// iter < 3 timestep = timestep * 1.3
// iter > 3 && iter < 7 timestep = timestep * 1.0

if (iter > 7)
{
    runTime += runTime.deltaT() * 0.7;
}
else if (iter < 3)
{
    runTime += runTime.deltaT() * 1.3;
}
else if (iter >= 3 && iter <= 7)
{
    runTime ++;
}