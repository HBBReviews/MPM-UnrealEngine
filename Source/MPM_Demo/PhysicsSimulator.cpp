// Fill out your copyright notice in the Description page of Project Settings.


#include "PhysicsSimulator.h"

#include "DrawDebugHelpers.h"

// Sets default values
APhysicsSimulator::APhysicsSimulator()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;
}

// Called when the game starts or when spawned
void APhysicsSimulator::BeginPlay()
{
	Super::BeginPlay();
	InitializePhysicsSimulator();
}

void APhysicsSimulator::InitializePhysicsSimulator()
{
	NumberOfCells = GridResolution * GridResolution;
	
	// Set up the current weights vectors.
	CurrentWeights.AddDefaulted(3);

	// Initialize all the points in a square.
	TArray<FVector2D> InitialParticlePositions = TArray<FVector2D>();
	
	const float SPACING = ParticleSpacing;
	const FVector2D BoxDimensions = FVector2D(BoxDimensionsSize);
	const FVector2D ScreenDimensions = FVector2D(GridResolution / 2.f);
	
	for (float X = BoxDimensions.X - BoxDimensions.X / 2; X < (ScreenDimensions.X + ScreenDimensions.X / 2); X += SPACING)
	{
		for (float Y = ScreenDimensions.Y - BoxDimensions.Y / 2; Y < ScreenDimensions.Y + BoxDimensions.Y / 2; Y += SPACING)
		{
			InitialParticlePositions.Add(FVector2D(X, Y));
		}
	}

	NumberOfParticles = InitialParticlePositions.Num();

	// Populate the array of particles and set the initial state.
	AllParticles.AddDefaulted(NumberOfParticles);

	for (int Index_Particle = 0; Index_Particle < NumberOfParticles; Index_Particle++)
	{
		FPhysicsParticle& NewParticle = AllParticles[Index_Particle];
		NewParticle.Position = InitialParticlePositions[Index_Particle];

		const float FirstRandom = FMath::RandRange(0.f, 1.f);
		const float SecondRandom = FMath::RandRange(0.f, 1.f);
		NewParticle.Velocity = FVector2D(FirstRandom - 0.5f, SecondRandom - 0.5f + 2.75f) * 0.5f;

		NewParticle.C = FPhysicsMatrix(0.0, 0.0, 0.0, 0.0);
		NewParticle.Mass = 1.f;
	}

	// Set up and populate the array of grid cells.
	CurrentGrid.AddDefaulted(NumberOfCells);
}

void APhysicsSimulator::SimulatePhysicsWithMPM(float DeltaTime)
{
	// Reset the grid's scratch pad.
	ResetGrid();

	RunParticleToGridSimulation();
	RunGridVelocity();
	RunGridToParticleSimulation();
}

// Called every frame
void APhysicsSimulator::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

	GEngine->AddOnScreenDebugMessage(-1, DeltaTime * 1.05f, FColor::Yellow, FString::Printf(TEXT("Current Force: %s"),
		*CurrentForceVector.ToString()));

	if ((ForceUpdateTimer += DeltaTime) >= TimeBetweenEachForceStrengthUpdate)
	{
		ForceUpdateTimer = 0.f;
		UpdateCurrentForce();
	}

	for (int Index_Iteration = 0; Index_Iteration < NumberOfIterationsPerTick; Index_Iteration++)
	{
		SimulatePhysicsWithMPM(DeltaTime);
	}

	Render();
}

void APhysicsSimulator::Render()
{
	// Render all of the particle locations.
	// ToDo: Convert this all to be Niagara.
	for (FPhysicsParticle ThisParticle : AllParticles)
	{
		DrawDebugPoint(GetWorld(), FVector(ThisParticle.Position, GetActorLocation().Z), DEBUG_PointSize,
			FColor::Blue, false, -1);
	}
}

void APhysicsSimulator::UpdateCurrentForce()
{
	const float RandomValue_1 = FMath::RandRange((float)-1.0, (float)1.0);
	const float RandomValue_2 = FMath::RandRange((float)-1.0, (float)1.0);

	CurrentForceVector = FVector2D(RandomValue_1 * ForceStrength, RandomValue_2 * ForceStrength);
}

void APhysicsSimulator::ResetGrid()
{
	for (int Index_Cell = 0; Index_Cell < NumberOfCells; ++Index_Cell)
	{
		FPhysicsGridCell& ThisGridCell = CurrentGrid[Index_Cell];
		ThisGridCell.Mass = 0.f;
		ThisGridCell.Velocity = FVector2D::ZeroVector;
	}
}

void APhysicsSimulator::RunParticleToGridSimulation()
{
	// Do the Particle to Grid logic.
	for (int Index_Particle = 0; Index_Particle < NumberOfParticles; Index_Particle++)
	{
		FPhysicsParticle& ThisParticle = AllParticles[Index_Particle];

		FVector2D CellIndex = FVector2D((int)ThisParticle.Position.X, (int)ThisParticle.Position.Y);
		const FVector2D CellDiff = ThisParticle.Position - CellIndex - FVector2D(0.5f);

		CurrentWeights[0] = FVector2D(0.5f * FMath::Pow(0.5f - CellDiff.X, 2), 0.5f * FMath::Pow(0.5f - CellDiff.Y, 2));
		CurrentWeights[1] = FVector2D(0.75f - FMath::Square(CellDiff.X), 0.75f - FMath::Square(CellDiff.Y));
		CurrentWeights[2] = FVector2D(0.5f * FMath::Square(0.5f + CellDiff.X), 0.5f * FMath::Square(0.5f + CellDiff.Y));

		// For all the surrounding 9 cells
		for (int Index_GridX = 0; Index_GridX < NUMBER_OF_GRID_CELLS_TO_SEARCH_IN_EACH_DIMENSION; Index_GridX++)
		{
			for (int Index_GridY = 0; Index_GridY < NUMBER_OF_GRID_CELLS_TO_SEARCH_IN_EACH_DIMENSION; Index_GridY++)
			{
				const float ThisWeight = CurrentWeights[Index_GridX].X * CurrentWeights[Index_GridX].Y;
				FVector2D CellLocation = FVector2D((int)(CellIndex.X + Index_GridX - 1), (int)(CellIndex.Y + Index_GridY - 1));

				const FVector2D CellDistance = CellLocation - ThisParticle.Position + FVector2D(0.5f);
				FVector2D Q = FVector2D(
					(ThisParticle.C.BackingArray[0] * CellDistance.X) + (ThisParticle.C.BackingArray[1] * CellDistance.X),
					(ThisParticle.C.BackingArray[2] * CellDistance.Y) + (ThisParticle.C.BackingArray[3] * CellDistance.Y)
					);

				float MassContribution = ThisWeight * ThisParticle.Mass;

				int CellIndex_Singular = (int)(CellLocation.X * GridResolution) + (int)CellLocation.Y;
				FPhysicsGridCell& ThisCell = CurrentGrid[CellIndex_Singular];

				ThisCell.Mass += MassContribution;
				ThisCell.Velocity += (MassContribution * (ThisParticle.Velocity + Q));
			}
		}
	}
}

void APhysicsSimulator::RunGridVelocity()
{
	for (int Index_Cell = 0; Index_Cell < NumberOfCells; Index_Cell++)
	{
		FPhysicsGridCell& ThisCell = CurrentGrid[Index_Cell];

		if (ThisCell.Mass > 0)
		{
			ThisCell.Velocity /= ThisCell.Mass;
			ThisCell.Velocity += FVector2D(0.f, GRAVITY_VALUE);

			const int CellIndex_X = (int)(Index_Cell / GridResolution);
			const int CellIndex_Y = (int)(Index_Cell % (int)GridResolution);

			if (CellIndex_X < 2 || CellIndex_X > GridResolution - 3)
			{
				ThisCell.Velocity.X = 0;
			}

			if (CellIndex_Y < 2 || CellIndex_Y > GridResolution - 3)
			{
				ThisCell.Velocity.Y = 0;
			}
		}
	}
}

void APhysicsSimulator::RunGridToParticleSimulation()
{
	for (int Index_Particle = 0; Index_Particle < AllParticles.Num(); Index_Particle++)
	{
		FPhysicsParticle& ThisParticle = AllParticles[Index_Particle];
		ThisParticle.Velocity = FVector2D::ZeroVector;

		FVector2D ParticleCellLocation = FVector2D((int)ThisParticle.Position.X, (int)ThisParticle.Position.Y);
		const FVector2D CellDifference = (ThisParticle.Position - ParticleCellLocation) - FVector2D(0.5f);

		CurrentWeights[0] = FVector2D(0.5f * FMath::Pow(0.5f - CellDifference.X, 2), 0.5f * FMath::Pow(0.5f - CellDifference.Y, 2));
		CurrentWeights[1] = FVector2D(0.75f - FMath::Square(CellDifference.X), 0.75f - FMath::Square(CellDifference.Y));
		CurrentWeights[2] = FVector2D(0.5f * FMath::Square(0.5f + CellDifference.X), 0.5f * FMath::Square(0.5f + CellDifference.Y));

		FPhysicsMatrix B = FPhysicsMatrix(0.f, 0.f, 0.f, 0.f);
		for (int GridX = 0; GridX < NUMBER_OF_GRID_CELLS_TO_SEARCH_IN_EACH_DIMENSION; ++GridX)
		{
			for (int GridY = 0; GridY < NUMBER_OF_GRID_CELLS_TO_SEARCH_IN_EACH_DIMENSION; ++GridY)
			{
				const float ThisWeight = CurrentWeights[GridX].X * CurrentWeights[GridX].Y;

				FVector2D ThisCellLocation = FVector2D((int)(ParticleCellLocation.X + GridX - 1),
					(int)(ParticleCellLocation.Y + GridY - 1));
				const FVector2D ThisCellDistance = ThisCellLocation - ThisParticle.Position + FVector2D(0.5f);

				const int CellIndex = (int)(ThisCellLocation.X * GridResolution) + (int)ThisCellLocation.Y;
				FPhysicsGridCell& ThisCell = CurrentGrid[CellIndex];

				FVector2D WeightedVelocity = ThisCell.Velocity * ThisWeight;

				// var term = math.float2x2(weighted_velocity * dist.x, weighted_velocity * dist.y);
				const FVector2D Term_A = WeightedVelocity * ThisCellDistance.X;
				const FVector2D Term_B = WeightedVelocity * ThisCellDistance.Y;

				// Manually set these up to be a 2x2 matrix.
				B.BackingArray[0] += Term_A.X;
				B.BackingArray[2] += Term_A.Y;
				
				B.BackingArray[1] += Term_B.X;
				B.BackingArray[3] += Term_B.Y;
				
				ThisParticle.Velocity += WeightedVelocity;
			}
		}

		ThisParticle.C = B;
		ThisParticle.C.BackingArray[0] *= 4;
		ThisParticle.C.BackingArray[1] *= 4;
		ThisParticle.C.BackingArray[2] *= 4;
		ThisParticle.C.BackingArray[3] *= 4;

		ThisParticle.Velocity += CurrentForceVector;

		// Update the particle with the current velocity.
		ThisParticle.Position += ThisParticle.Velocity;
		
		ThisParticle.Position.X = FMath::Clamp((int)ThisParticle.Position.X, 1, (int)GridResolution - 2);
		ThisParticle.Position.Y = FMath::Clamp((int)ThisParticle.Position.Y, 1, (int)GridResolution - 2);

		// // constructing affine per-particle momentum matrix from APIC / MLS-MPM.
		// // see APIC paper (https://web.archive.org/web/20190427165435/https://www.math.ucla.edu/~jteran/papers/JSSTS15.pdf), page 6
		// // below equation 11 for clarification. this is calculating C = B * (D^-1) for APIC equation 8,
		// // where B is calculated in the inner loop at (D^-1) = 4 is a constant when using quadratic interpolation functions
		// float2x2 B = 0;
		// for (uint gx = 0; gx < 3; ++gx) {
		// 	for (uint gy = 0; gy < 3; ++gy) {
		// 		float weight = weights[gx].x * weights[gy].y;
		//
		// 		uint2 cell_x = math.uint2(cell_idx.x + gx - 1, cell_idx.y + gy - 1);
		// 		int cell_index = (int)cell_x.x * grid_res + (int)cell_x.y;
		//
		// 		float2 dist = (cell_x - p.x) + 0.5f;
		// 		float2 weighted_velocity = grid[cell_index].v * weight;
		//
		// 		// APIC paper equation 10, constructing inner term for B
		// 		var term = math.float2x2(weighted_velocity * dist.x, weighted_velocity * dist.y);
		//
		// 		B += term;
		//
		// 		p.v += weighted_velocity;
		// 	}
		// }
		// p.C = B * 4;
		//
		// // advect particles
		// p.x += p.v * dt;
		//
		// // safety clamp to ensure particles don't exit simulation domain
		// p.x = math.clamp(p.x, 1, grid_res - 2);
	}
}
