void AllocArrays ();
void Init_F ();
void CBD_Step ();
void SetupJob ();
void SingleStep ();
void PrintSummary ();
void ReadInput ();
void writeProps();
void feneForce ();
void InitCoords ();
void computeForces (int stage);
void PrintElapsedTime(chrono::steady_clock::time_point start);
void writePropsHeader();
void averageProps();
void cellDivision();
void cellDeath();
void wall();
