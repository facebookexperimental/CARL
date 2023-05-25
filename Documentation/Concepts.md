# CARL Conceptual Overview

CARL is not a particularly complicated technology, but the ideas that 
underpin it can be specialized and worth describing at a conceptual
level. This document will begin with a brief discussion of the objectives 
and perspectives informing CARL's approach to action recognition, then 
discuss the concepts which emerge from this approach and are directly reflected in the implementation.

## CARL's Objectives and Perspectives

In purpose and style, CARL is perhaps most closely related to the 
[$-family of recognizers](https://depts.washington.edu/acelab/proj/dollar/impact.html).
This is a family of algorithms and approaches which characteristically
tackle various manifestations of the action recognition problem using
"classical" techniques to directly characterize and compare actions without
reliance on large datasets or advanced statistical models, thus allowing
usage in a much more flexible set of circumstances. Similar sensibilities 
also inform CARL's approach, which is designed to be as accessible as 
possible and to require as little effort to operate as possible.

The most fundamental goal of CARL has always been to allow actions to be
defined by simply _doing_ them. This is how people (and even animals) will
naturally define actions when teaching them to one another -- "Here, just
watch what I'm doing" -- and thus it seems appropriate to express action
definitions for machines in a similar way. This concept of action 
definition -- that an action is best defined by examples of it being done
correctly -- also means that new definitions can be created by anyone
who can perform the desired motion, without requiring substantial tooling,
large datasets, or specialized expertise. This opens the door for 
learning/updating action definitions "on the fly," which has significant
implications for tailoring action recognition to individual users.

Given the above perspective on action recognition, the core concepts of CARL naturally emerge.
- An action is defined by a **definition**.
- A **definition** consists of one or more **example**s of the action being performed correctly.
- An **example** consists of a **recording** of user state during a timeframe when the user was performing an action.
- A **recording** is a time series of **input sample**s.
- An **input sample** is an instanteous "snapshot" of the relevant state from which actions can be defined and recognized.

These concepts, plus the **session**s and **recognizer**s which are needed 
to perform recognition, are described in more detail in the following 
sections.

## Input Sample

An `InputSample` is CARL's atomic unit of input state: it is a "snapshot" of
every piece of relevant information available at a certain point in time. 
For example, for an `InputSample` describing an XR participant wearing an
HMD and using hand tracking, an `InputSample` might contain the 3D pose 
(position and orientation) of the HMD itself and every joint in each tracked
hand at a single point in time. Since not all of this data may be available 
all the time (for example, a hand might not be tracked at times), separable
components of this data are optional. The only never-optional element of an
`InputSample` is its timestamp, which is expressed in seconds.

Much of the data contained in an `InputSample` consists of 3D poses. By 
convention, all described poses should be expressed in the same 
left-handed Y-up coordinate space. Other conventions are not yet supported
and behavior on data with such conventions is undefined.

## Recording

A `Recording` is a continuous sequence of `InputSample`s, ordered from 
earliest to latest according to their timestamps. 

Note that the timestamps of `InputSample`s in a `Recording` do not have to 
be relative to any particular reference time (epoch time, etc.); all that 
is necessary for a `Recording` to be valid is that the timestamps of all 
`InputSample`s be relative to the _same_ reference time. Also note that 
`Recording`s do nothave any specific framerate and are not required to 
contain `InputSample`sat any particular (or even consistent) frequency.

## Example

An `Example` is a `Recording` which is known to contain a manifestation of
a semantically-definite action. Specifically, the action is known to have 
"happened" between a given `startTimestamp` and a given `endTimestamp`, 
such that `InputSample`s in the `Recording` between those two timestamps 
constitute an "example" of that action occuring.

Note that `Example` does not reference any form of descriptor or 
definition; in fact, `Example` contains no explicit reference to any 
particular semantically-definite action at all. Colloquially speaking,
an `Example` specifies nothing about _what_ happened, but merely an 
example of _something_ happening.

## Definition

A `Definition` is a set of `Example`s which are all known to contain 
manifestations of the same semantically-definite action. It also identifies
a `DescriptorType` specifying what kind of descriptor should be used to 
evaluate the defined action. In this way, `Definition` serves as CARL's 
fundamental concretization of the concept of semantically-definite actions.

## Session

A `Session` can be thought of as a representation of a time period when 
actions can be recognized. Contractually, the `Session` itself contains
nothing and exposes nothing; it simply functions as a receptacle into which
`InputSample`s should be placed, and as a resource required to build
`Recognizer`s.

## Recognizer

A `Recognizer` is a construct which observes a `Session` over time (i.e.,
as new `InputSample`s are added) and expresses a "confidence" that a
semantically-definite action defined by a `Definition` has been recognized.
This "confidence" is a floating-point number (currently a `double`) which
is large (near or above 1) when the action in question is recognized and 
small (near or below 0) when it is not recognized. For simple scenarios, 
this number can simply be clamped to the [0,1] range and treated like a
probability (though it is not probabilistic in nature). More sophsticated 
usages like automatic example segmentation currently leverage unbounded
confidence scoring; however, this is an evolving usage that will likely
change/simplify in the near future.
