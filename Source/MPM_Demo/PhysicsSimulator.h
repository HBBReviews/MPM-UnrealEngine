// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "PhysicsSimulator.generated.h"

USTRUCT()
struct FPhysicsMatrix
{
	GENERATED_BODY()

	UPROPERTY(VisibleInstanceOnly)
	TArray<float> BackingArray;

	UPROPERTY(VisibleInstanceOnly)
	FVector2D Size = FVector2D::ZeroVector;
	
	FPhysicsMatrix(float Value1, float Value2, float Value3, float Value4)
	{
		Size = FVector2D(2, 2);
		BackingArray.AddDefaulted(4);

		BackingArray[0] = Value1;
		BackingArray[1] = Value2;
		BackingArray[2] = Value3;
		BackingArray[3] = Value4;
	}

	FPhysicsMatrix()
	{
		FPhysicsMatrix(0.f, 0.f, 0.f, 0.f);
	}
};

USTRUCT()
struct FPhysicsParticle
{
	GENERATED_BODY()

	UPROPERTY(VisibleInstanceOnly)
	FVector2D Position = FVector2D::ZeroVector;

	UPROPERTY(VisibleInstanceOnly)
	FVector2D Velocity = FVector2D::ZeroVector;

	// Interpret this as a 2x2 Matrix.
	UPROPERTY(VisibleInstanceOnly)
	FPhysicsMatrix C = FPhysicsMatrix();

	UPROPERTY(VisibleInstanceOnly)
	float Mass = 0.f;

	UPROPERTY(VisibleInstanceOnly)
	float Padding = 0.f;
};

USTRUCT()
struct FPhysicsGridCell
{
	GENERATED_BODY()

	UPROPERTY()
	FVector2D Velocity = FVector2D::ZeroVector;

	UPROPERTY()
	float Mass = 0.f;

	UPROPERTY()
	float Padding = 0.f;
};

UCLASS()
class MPM_DEMO_API APhysicsSimulator : public AActor
{
	GENERATED_BODY()

	static constexpr float GRAVITY_VALUE = -0.05f;
	static constexpr int NUMBER_OF_GRID_CELLS_TO_SEARCH_IN_EACH_DIMENSION = 3;

public:	
	// Sets default values for this actor's properties
	APhysicsSimulator();

	UPROPERTY(EditDefaultsOnly, Category = "Physics Simulator|Settings")
	int NumberOfIterationsPerTick = 1;

	UPROPERTY(EditDefaultsOnly, Category = "Physics Simulator|Settings")
	float ParticleSpacing = 1.0;

	UPROPERTY(EditDefaultsOnly, Category = "Physics Simulator|Settings")
	int BoxDimensionsSize = 16;

	UPROPERTY(EditDefaultsOnly, Category = "Physics Simulator|Settings")
	float GridResolution = 64;

	UPROPERTY(EditDefaultsOnly, Category = "Physics Simulator|Settings")
	float ForceStrength = 5.f;

	UPROPERTY(EditDefaultsOnly, Category = "Physics Simulator|Settings")
	float DEBUG_PointSize = 30.f;

	UPROPERTY(EditDefaultsOnly, Category = "Physics Simulator|Settings")
	float TimeBetweenEachForceStrengthUpdate = 5.f;

	UPROPERTY(VisibleInstanceOnly, Category = "Physics Simulator|Particles")
	int NumberOfParticles = 0;

	UPROPERTY(VisibleInstanceOnly, Category = "Physics Simulator|Particles")
	TArray<FPhysicsParticle> AllParticles;

	UPROPERTY(VisibleInstanceOnly, Category = "Physics Simulator|Grid")
	TArray<FPhysicsGridCell> CurrentGrid;

	UPROPERTY(VisibleInstanceOnly, Category = "Physics Simulator|Grid")
	int NumberOfCells = 0;
	
	UPROPERTY()
	TArray<FVector2D> CurrentWeights;

	UPROPERTY()
	FVector2D CurrentForceVector = FVector2D::ZeroVector;
	
	float ForceUpdateTimer = 0.f;

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

	void InitializePhysicsSimulator();
	void SimulatePhysicsWithMPM(float DeltaTime);
	void Render();

	void UpdateCurrentForce();

	void ResetGrid();
	
	void RunParticleToGridSimulation();
	void RunGridVelocity();
	void RunGridToParticleSimulation();
	
public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

};
